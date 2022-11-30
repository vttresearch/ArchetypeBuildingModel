#=
    create_spineopt_input.jl

Functions and structs for handling input data for SpineOpt based on the archetype building model.
=#

"""
    SpineOptInput

Create and store the input data for the SpineOpt energy system model.

Contains the following fields:
- `node::ObjectClass`: Contains all the `node`s in the building models, created based on [`AbstractNode`](@ref)s.
- `unit::ObjectClass`: Contains all the `unit`s in the building models, created based on [`AbstractProcess`](@ref)es.
- `node__node::RelationshipClass`: Defines the heat transfer coefficients between the `node`s.
- `unit__from_node::RelationshipClass`: Defines `unit` input flow properties.
- `unit__to_node::RelationshipClass`: Defines `unit` output flow properties.
- `unit__node__node::RelationshipClass`: Defines `unit` conversion properties.
"""
struct SpineOptInput <: ModelInput
    node::ObjectClass
    unit::ObjectClass
    node__node::RelationshipClass
    unit__from_node::RelationshipClass
    unit__to_node::RelationshipClass
    unit__node__node::RelationshipClass
    function SpineOptInput()
        node = ObjectClass(:node, Array{ObjectLike,1}())
        unit = ObjectClass(:unit, Array{ObjectLike,1}())
        node__node =
            RelationshipClass(:node__node, [:node, :node], Array{RelationshipLike,1}())
        unit__from_node =
            RelationshipClass(:unit__from_node, [:unit, :node], Array{RelationshipLike,1}())
        unit__to_node =
            RelationshipClass(:unit__to_node, [:unit, :node], Array{RelationshipLike,1}())
        unit__node__node = RelationshipClass(
            :unit__node__node,
            [:unit, :node, :node],
            Array{RelationshipLike,1}(),
        )
        new(node, unit, node__node, unit__from_node, unit__to_node, unit__node__node)
    end
end


"""
    SpineOptInput(
        results::Dict{Object,ArchetypeBuildingResults};
        mod::Module = @__MODULE__,
    )

Create [`SpineOptInput`](@ref) based on given archetype building results.

NOTE! The `mod` keyword changes from which Module data is accessed from,
`@__MODULE__` by default.

Essentially, performs the following steps:
1. Initialize an empty [`SpineOptInput`](@ref).
2. Loop over the given `results`, and [`add_archetype_to_input!`](@ref) one by one.
"""
function SpineOptInput(
    results::Dict{Object,ArchetypeBuildingResults};
    mod::Module = @__MODULE__,
)
    spineopt = SpineOptInput()
    for result in values(results)
        add_archetype_to_input!(spineopt, result; mod = mod)
    end
    return spineopt
end


"""
    add_archetype_to_input!(
        spineopt::SpineOptInput,
        result::ArchetypeBuildingResults;
        mod::Module = @__MODULE__,
    )

Process and add the desired archetype building `result` into the `spineopt` input.

NOTE! The `mod` keyword changes from which Module data is accessed from,
`@__MODULE__` by default.

This is a very long and rather complicated function,
which essentially translates the information contained in the `result` [`ArchetypeBuildingResults`](@ref)
into the data structure in `spineopt` [`SpineOptInput`](@ref),
so that it is understood by the SpineOpt energy system model.
The key steps taken by this function are summarized below:
1. Map the ArchetypeBuildingModel `building_node` objects to unique SpineOpt `node` objects.
2. Identify the necessary *system link nodes*, and add them into the set of SpineOpt `node` objects with the desired names.
3. Map the ArchetypeBuildingModel `building_process` objects to unique SpineOpt `unit` objects.
4. Determine SpineOpt `node` parameters based on the [`AbstractNode`](@ref) properties.
5. Determine SpineOpt `node__node` parameters based on the [`AbstractNode`](@ref) heat transfer coefficients.
6. Determine SpineOpt `unit__from_node` and `unit__to_node` parameters based on the [`AbstractProcess`](@ref) maximum flow parameters.
7. Determine SpineOpt `unit__node__node` parameters based on the [`AbstractProcess`](@ref) properties.
"""
function add_archetype_to_input!(
    spineopt::SpineOptInput,
    result::ArchetypeBuildingResults;
    mod::Module = @__MODULE__,
)
    # Map `building_node` objects to unique `node` objects.
    n_map = Dict(
        node => Object(
            Symbol(string(result.archetype.archetype.name) * "__" * string(node.name)),
            :node,
        ) for node in keys(result.archetype.abstract_nodes)
    )
    # Identify the relevant system link nodes and map them to new nodes with the desired names.
    sys_link_nodes = mod.building_archetype__system_link_node(
        building_archetype = result.archetype.archetype,
    )
    # This is a bit complicated, as we have to avoid creating identical objects.
    merge!(
        n_map,
        Dict(
            sys_link_node =>
                isnothing(
                    spineopt.node(
                        mod.node_name(
                            building_archetype = result.archetype.archetype,
                            building_node = sys_link_node,
                        ),
                    ),
                ) ?
                Object(
                    mod.node_name(
                        building_archetype = result.archetype.archetype,
                        building_node = sys_link_node,
                    ),
                    :node,
                ) :
                spineopt.node(
                    mod.node_name(
                        building_archetype = result.archetype.archetype,
                        building_node = sys_link_node,
                    ),
                ) for sys_link_node in sys_link_nodes
        ),
    )
    # Map `building_process` objects to unique `unit` objects.
    u_map = Dict(
        process => Object(
            Symbol(string(result.archetype.archetype.name) * "__" * string(process.name)),
            :unit,
        ) for process in keys(result.archetype.abstract_processes)
    )

    # Add the necessary `unit` objects
    add_objects!(spineopt.unit, collect(values(u_map)))

    # Add the necessary `node` objects and their parameters
    node_param_dict = Dict(
        n_map[node] => Dict(
            :demand => parameter_value(-1 * abs_node.external_load),
            :fix_node_state => parameter_value(
                TimeSeries(
                    [first(result.temperatures[node].indexes)],
                    [result.initial_temperatures[node]],
                    false,
                    false,
                ),
            ),
            :frac_state_loss =>
                parameter_value(abs_node.self_discharge_coefficient_W_K),
            :has_state => parameter_value(true),
            :node_state_cap => parameter_value(abs_node.maximum_temperature_K),
            :node_state_min => parameter_value(abs_node.minimum_temperature_K),
            :state_coeff => parameter_value(abs_node.thermal_mass_Wh_K),
        ) for (node, abs_node) in result.archetype.abstract_nodes
    )
    add_object_parameter_values!(spineopt.node, node_param_dict)
    merge!(
        spineopt.node.parameter_defaults,
        Dict(
            :demand => parameter_value(0.0),
            :fix_node_state => parameter_value(nothing),
            :frac_state_loss => parameter_value(0.0),
            :has_state => parameter_value(false),
            :node_state_cap => parameter_value(nothing),
            :node_state_min => parameter_value(0.0),
            :state_coeff => parameter_value(1.0),
        ),
    )

    # Add the necessary relationships and their parameters
    # `node__node` based on heat transfer coefficients
    nn_param_dict = Dict(
        (n_map[n1], n_map[n2]) => Dict(:diff_coeff => parameter_value(v)) for
        (n1, abs_n1) in result.archetype.abstract_nodes for
        (n2, v) in abs_n1.heat_transfer_coefficients_W_K
    )
    add_relationships!(
        spineopt.node__node,
        [(node1 = n1, node2 = n2) for (n1, n2) in keys(nn_param_dict)],
    )
    merge!(spineopt.node__node.parameter_values, nn_param_dict)
    spineopt.node__node.parameter_defaults[:diff_coeff] = parameter_value(0.0)

    # `unit__from_node` based on maximum flow parameters.
    ufn_param_dict = Dict(
        (u_map[p], n_map[n]) => Dict(:unit_capacity => parameter_value(abs(v))) for
        (p, abs_p) in result.archetype.abstract_processes for
        ((d, n), v) in abs_p.maximum_flows_W if (
            v >= 0 && d == mod.direction(:from_node) ||
            v < 0 && d == mod.direction(:to_node)
        )
    )
    add_relationships!(
        spineopt.unit__from_node,
        [(unit = u, from_node = n) for (u, n) in keys(ufn_param_dict)],
    )
    merge!(spineopt.unit__from_node.parameter_values, ufn_param_dict)
    spineopt.unit__from_node.parameter_defaults[:unit_capacity] = parameter_value(nothing)

    # `unit__to_node` based on maximum flow parameters.
    utn_param_dict = Dict(
        (u_map[p], n_map[n]) => Dict(:unit_capacity => parameter_value(abs(v))) for
        (p, abs_p) in result.archetype.abstract_processes for
        ((d, n), v) in abs_p.maximum_flows_W if (
            v >= 0 && d == mod.direction(:to_node) ||
            v < 0 && d == mod.direction(:from_node)
        )
    )
    add_relationships!(
        spineopt.unit__to_node,
        [(unit = u, to_node = n) for (u, n) in keys(utn_param_dict)],
    )
    merge!(spineopt.unit__to_node.parameter_values, utn_param_dict)
    spineopt.unit__to_node.parameter_defaults[:unit_capacity] = parameter_value(nothing)

    # `unit__node__node` based on maximum flow parameters as well.
    unn_param_dict = Dict(
        (u_map[p], n_map[n1], n_map[n2]) => Dict(
            :fix_ratio_in_out_unit_flow => parameter_value(
                v2 >= 0 ?
                1 / (
                    abs_p.coefficient_of_performance isa Union{TimeSeries,TimePattern} ?
                    timedata_operation(abs, abs_p.coefficient_of_performance) :
                    abs(abs_p.coefficient_of_performance)
                ) : nothing,
            ),
            :fix_ratio_in_in_unit_flow => parameter_value(
                v2 < 0 ?
                1 / (
                    abs_p.coefficient_of_performance isa Union{TimeSeries,TimePattern} ?
                    timedata_operation(abs, abs_p.coefficient_of_performance) :
                    abs(abs_p.coefficient_of_performance)
                ) : nothing,
            ),
        ) for (p, abs_p) in result.archetype.abstract_processes for
        ((d1, n1), v1) in abs_p.maximum_flows_W for
        ((d2, n2), v2) in abs_p.maximum_flows_W if
        (d1 == mod.direction(:from_node) && d2 == mod.direction(:to_node) && n1 != n2)
    )
    add_relationships!(
        spineopt.unit__node__node,
        [(unit = u, node1 = n1, node2 = n2) for (u, n1, n2) in keys(unn_param_dict)],
    )
    merge!(spineopt.unit__node__node.parameter_values, unn_param_dict)
    spineopt.unit__node__node.parameter_defaults[:fix_ratio_in_out_unit_flow] =
        parameter_value(nothing)
    spineopt.unit__node__node.parameter_defaults[:fix_ratio_in_in_unit_flow] =
        parameter_value(nothing)
end
