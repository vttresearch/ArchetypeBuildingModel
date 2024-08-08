#=
    create_spineopt_input.jl

Functions and structs for handling input data for SpineOpt based on the archetype building model.
=#

"""
    SpineOptInput(url::Union{String,Dict})

Link to the input data store for the SpineOpt energy system model.

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
    function SpineOptInput(url::Union{String,Dict}) # Fetch and link the database structure from url.
        m = Module() # Create a separate module to load SpineOpt data store structure into.
        using_spinedb(url, m)
        new(
            m.node,
            m.unit,
            m.node__node,
            m.unit__from_node,
            m.unit__to_node,
            m.unit__node__node
        )
    end
end


"""
    SpineOptInput(
        url::Union{String,Dict},
        results::Dict{Object,ArchetypeBuildingResults};
        mod::Module = @__MODULE__,
    )

Create [`SpineOptInput`](@ref) based on given archetype building results.

NOTE! The `mod` keyword changes from which Module data is accessed from,
`@__MODULE__` by default.

Essentially, performs the following steps:
1. Link to [`SpineOptInput`](@ref) at `url`.
2. Loop over the given `results`, and [`add_archetype_to_input!`](@ref) one by one.
"""
function SpineOptInput(
    url::Union{String,Dict},
    results::Dict{Object,ArchetypeBuildingResults};
    mod::Module=@__MODULE__
)
    spineopt = SpineOptInput(url)
    for result in values(results)
        add_archetype_to_input!(spineopt, result; mod=mod)
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
    mod::Module=@__MODULE__
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
        building_archetype=result.archetype.archetype,
    )
    # This is a bit complicated, as we have to avoid creating identical objects.
    merge!(
        n_map,
        Dict(
            sys_link_node =>
                isnothing(
                    spineopt.node(
                        mod.node_name(
                            building_archetype=result.archetype.archetype,
                            building_node=sys_link_node,
                        ),
                    ),
                ) ?
                Object(
                    mod.node_name(
                        building_archetype=result.archetype.archetype,
                        building_node=sys_link_node,
                    ),
                    :node,
                ) :
                spineopt.node(
                    mod.node_name(
                        building_archetype=result.archetype.archetype,
                        building_node=sys_link_node,
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
            :demand => parameter_value(-1 * abs_node.external_load_kW),
            :initial_node_state => parameter_value(result.initial_temperatures[node]),
            :frac_state_loss =>
                parameter_value(abs_node.self_discharge_coefficient_kW_K),
            :has_state => parameter_value(true),
            :node_state_cap => parameter_value(abs_node.maximum_temperature_K),
            :node_state_min => parameter_value(abs_node.minimum_temperature_K),
            :state_coeff => parameter_value(abs_node.thermal_mass_kWh_K),
        ) for (node, abs_node) in result.archetype.abstract_nodes
    )
    add_object_parameter_values!(spineopt.node, node_param_dict; merge_values=true)

    # Add the necessary relationships and their parameters
    # `node__node` based on heat transfer coefficients
    nn_param_dict = Dict(
        (n_map[n1], n_map[n2]) => Dict(:diff_coeff => parameter_value(v)) for
        (n1, abs_n1) in result.archetype.abstract_nodes for
        (n2, v) in abs_n1.heat_transfer_coefficients_kW_K
    )
    add_relationship_parameter_values!(spineopt.node__node, nn_param_dict; merge_values=true)

    # `unit__from_node` based on maximum flow parameters.
    ufn_param_dict = Dict(
        (u_map[p], n_map[n]) => Dict(:unit_capacity => parameter_value(abs(v))) for
        (p, abs_p) in result.archetype.abstract_processes for
        ((d, n), v) in abs_p.maximum_flows if (
            v >= 0 && d == mod.direction(:from_node) ||
            v < 0 && d == mod.direction(:to_node)
        )
    )
    add_relationship_parameter_values!(spineopt.unit__from_node, ufn_param_dict; merge_values=true)

    # `unit__to_node` based on maximum flow parameters.
    utn_param_dict = Dict(
        (u_map[p], n_map[n]) => Dict(:unit_capacity => parameter_value(abs(v))) for
        (p, abs_p) in result.archetype.abstract_processes for
        ((d, n), v) in abs_p.maximum_flows if (
            v >= 0 && d == mod.direction(:to_node) ||
            v < 0 && d == mod.direction(:from_node)
        )
    )
    add_relationship_parameter_values!(spineopt.unit__to_node, utn_param_dict; merge_values=true)

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
        ((d1, n1), v1) in abs_p.maximum_flows for
        ((d2, n2), v2) in abs_p.maximum_flows if
        (d1 == mod.direction(:from_node) && d2 == mod.direction(:to_node) && n1 != n2)
    )
    add_relationship_parameter_values!(spineopt.unit__node__node, unn_param_dict; merge_values=true)
end
