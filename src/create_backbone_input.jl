#=
    create_backbone_input.jl

Functions for creating Backbone input data from the archetype buildings.
=#

"""
    BackboneInput(url::{String,Dict})

Link to the input data store for the Backbone energy system model.

Contains the following fields:
- `boundary::ObjectClass`: Contains the `upwardLimit` and `downwardLimit` settings.
- `effLevel::ObjectClass`: All possible efficiency representation levels in Backbone.
- `effSelector::ObjectClass`: Contains the `directOff` efficiency representation used by the building model.
- `grid::ObjectClass`: Contains the `building` `grid` for building `node`s.
- `io::ObjectClass`: Contains the `input` and `output` indicators for `unit`s.
- `node::ObjectClass`: Contains all the `node`s in the building models, created based on [`AbstractNode`](@ref)s.
- `unit::ObjectClass`: Contains all the `unit`s in the building models, created based on [`AbstractProcess`](@ref)es.
- `unittype::ObjectClass`: Contains a `HVAC` type for all building `unit`s.
- `effLevel__effSelector__unit::RelationshipClass`: Attributes the `directOff` efficiency representation for all `unit`s for all efficiency representation levels.
- `grid__node::RelationshipClass`: Connects all `node`s into the `building` `grid`.
- `grid__node__boundary::RelationshipClass`: Contains the `upwardLimit` and `downwardLimit` parameters for all `node`s.
- `grid__node__node::RelationshipClass`: Contains the `diffCoeff` diffusion parameters between the `node`s.
- `grid__node__unit__io::RelationshipClass`: Defines how the `unit`s interact with the `node`s, and contains the necessary parameters.
- `unit__unittype::RelationshipClass`: Indicates all `unit`s as of type `HVAC`
"""
struct BackboneInput <: ModelInput
    boundary::ObjectClass
    effLevel::ObjectClass
    effSelector::ObjectClass
    grid::ObjectClass
    io::ObjectClass
    node::ObjectClass
    unit::ObjectClass
    unittype::ObjectClass
    effLevel__effSelector__unit::RelationshipClass
    grid__node::RelationshipClass
    grid__node__boundary::RelationshipClass
    grid__node__node::RelationshipClass
    grid__node__unit__io::RelationshipClass
    unit__unittype::RelationshipClass
    function BackboneInput(url::Union{String,Dict}) # Fetch and link the database structure from url.
        m = Module() # Create a separate module to load Backbone data store structure into.
        using_spinedb(url, m)
        new(
            m.boundary,
            m.effLevel,
            m.effSelector,
            m.grid,
            m.io,
            m.node,
            m.unit,
            m.unittype,
            m.effLevel__effSelector__unit,
            m.grid__node,
            m.grid__node__boundary,
            m.grid__node__node,
            m.grid__node__unit__io,
            m.unit__unittype,
        )
    end
end


"""
    BackboneInput(
        url::Union{String,Dict},
        results::Dict{Object,ArchetypeBuildingResults};
        mod::Module = @__MODULE__
    )

Create [`BackboneInput`](@ref) based on given archetype building results.

NOTE! The `mod` keyword changes from which Module data is accessed from,
`@__MODULE__` by default.

Essentially, performs the following steps:
1. Link to [`BackboneInput`](@ref) at `url`.
2. Loop over the given `results`, and [`add_archetype_to_input!`](@ref) one by one.
3. Calculate and [`add_system_link_node_parameters!`](@ref).
"""
function BackboneInput(
    url::Union{String,Dict},
    results::Dict{Object,ArchetypeBuildingResults};
    mod::Module=@__MODULE__
)
    backbone = BackboneInput(url)
    for result in values(results)
        add_archetype_to_input!(backbone, result; mod=mod)
    end
    add_system_link_node_parameters!(backbone, results; mod=mod)
    return backbone
end


"""
    add_archetype_to_input!(
        backbone::BackboneInput,
        result::ArchetypeBuildingResults;
        mod::Module = @__MODULE__,
    )

Process and add the desired archetype building `result` into the `backbone` input.

NOTE! The `mod` keyword changes from which Module data is accessed from,
`@__MODULE__` by default.

This is a very long and rather complicated function,
which essentially translates the information contained in the `result` [`ArchetypeBuildingResults`](@ref)
into the data structure in `backbone` [`BackboneInput`](@ref),
so that it is understood by the Backbone energy system model.
The key steps taken by this function are summarized below:
1. Map the ArchetypeBuildingModel `direction` objects to Backbone `io` objects.
2. Create a `grid` object representing the `archetype` being processed and map the grids.
3. Map the ArchetypeBuildingModel `building_node` objects to unique Backbone `node` objects.
4. Identify the necessary *system link nodes*, and add them into the set of Backbone `node` objects with the desired names.
5. Include all building `node`s to the archetype `grid`, and connect the *system link nodes* to their desired `grid`s.
6. Map the ArchetypeBuildingModel `building_process` objects to unique Backbone `unit` objects.
7. Determine Backbone `unit` parameters based on [`AbstractProcess`](@ref) properties.
8. Set the `directOff` efficiency representation for all `unit`s.
9. Determine Backbone `grid__node` parameters based on [`AbstractNode`](@ref) properties.
10. Determine Backbone `grid__node__boundary` parameters based on [`AbstractNode`](@ref) maximum and minimum permitted temperatures.
11. Determine Backbone `grid__node__node` parameters based on [`AbstractNode`](@ref) heat transfer coefficients.
12. Determine Backbone `grid__node__unit__io` parameters based on [`AbstractProcess`](@ref) maximum flows.
13. Set the `HVAC` `unittype` for every `unit`.
"""
function add_archetype_to_input!(
    backbone::BackboneInput,
    result::ArchetypeBuildingResults;
    mod::Module=@__MODULE__
)
    # Map `direction` objects to `io` objects.
    io_map = Dict(
        mod.direction(:to_node) => backbone.io(:output),
        mod.direction(:from_node) => backbone.io(:input),
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
                    backbone.node(
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
                backbone.node(
                    mod.node_name(
                        building_archetype=result.archetype.archetype,
                        building_node=sys_link_node,
                    ),
                ) for sys_link_node in sys_link_nodes
        ),
    )
    # Map `building_node` objects to their corresponding `grid
    arch_grid = Object(result.archetype.archetype.name, :grid)
    g_map = Dict(node => arch_grid for node in keys(result.archetype.abstract_nodes))
    # Map the system link nodes to their desired grids.
    # This is a bit complicated, as we have to avoid creating identical objects.
    merge!(
        g_map,
        Dict(
            sys_link_node =>
                isnothing(
                    backbone.grid(
                        mod.grid_name(
                            building_archetype=result.archetype.archetype,
                            building_node=sys_link_node,
                        ),
                    ),
                ) ?
                Object(
                    mod.grid_name(
                        building_archetype=result.archetype.archetype,
                        building_node=sys_link_node,
                    ),
                    :grid,
                ) :
                backbone.grid(
                    mod.grid_name(
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
    # Add various auxiliary `grid` parameters for results processing convenience
    g_param_dict = Dict(
        arch_grid => Dict(
            :ambient_temperature_K => parameter_value(
                timeseries_to_backbone_map(
                    result.archetype.weather_data.ambient_temperature_K,
                ),
            ),
            :ground_temperature_K => parameter_value(
                timeseries_to_backbone_map(
                    result.archetype.weather_data.ground_temperature_K,
                ),
            ),
            :number_of_buildings =>
                parameter_value(result.archetype.scope_data.number_of_buildings),
            :average_gross_floor_area_m2_per_building => parameter_value(
                result.archetype.scope_data.average_gross_floor_area_m2_per_building,
            ),
            Symbol("archetype.", result.archetype.archetype.name) => parameter_value(1),
            Symbol("scope.", result.archetype.scope.name) => parameter_value(1),
            Symbol("fabrics.", result.archetype.fabrics.name) => parameter_value(1),
            Symbol("systems.", result.archetype.systems.name) => parameter_value(1),
            Symbol("loads.", result.archetype.loads.name) => parameter_value(1),
            Symbol("weather.", result.archetype.weather.name) => parameter_value(1),
        ),
    )
    add_object_parameter_values!(backbone.grid, g_param_dict; merge_values=true)

    # Create the necessary objects and their parameters
    # Add the `grid` and `node` objects.
    add_objects!(backbone.grid, unique(values(g_map)))
    add_objects!(backbone.node, unique(values(n_map)))

    # Add `unit` objects and its parameters.
    u_param_dict = Dict(
        u_map[p] => Dict(
            :availability => parameter_value(1.0),
            :efficiency => parameter_value( # This is only used if COP isn't a time series
                abs_p.coefficient_of_performance isa Union{TimeSeries,TimePattern} ? nothing :
                Map([:op00, :eff00], [1.0, abs(abs_p.coefficient_of_performance)]),
            ),
            :efficiency_ts => parameter_value( # This is used if COP is a time series
                abs_p.coefficient_of_performance isa Union{TimeSeries,TimePattern} ?
                Map(
                    [:eff00],
                    [timeseries_to_backbone_map(abs_p.coefficient_of_performance)],
                ) : nothing,
            ),
            :unitCount => parameter_value(abs_p.number_of_processes),
            :useTimeseries => parameter_value( # Set flag for time series
                abs_p.coefficient_of_performance isa Union{TimeSeries,TimePattern},
            ),
        ) for (p, abs_p) in result.archetype.abstract_processes
    )
    add_object_parameter_values!(backbone.unit, u_param_dict; merge_values=true)

    # Create and add the necessary relationships
    # `effLevel__effSelector__unit` assumed static `directOff`.
    add_relationships!(
        backbone.effLevel__effSelector__unit,
        [
            (effLevel=lvl, effSelector=backbone.effSelector(:directOff), unit=u) for
            lvl in backbone.effLevel() for u in backbone.unit()
        ],
    )

    # `grid__node` and its parameters from `abstract_nodes`
    gn_param_dict = Dict(
        (grid=g_map[n], node=n_map[n]) => Dict(
            :boundStart => parameter_value(true),
            :energyStoredPerUnitOfState => parameter_value(abs_n.thermal_mass_kWh_K),
            :nodeBalance => parameter_value(true),
            :influx =>
                parameter_value(timeseries_to_backbone_map(abs_n.external_load_kW)),
            :selfDischargeLoss => parameter_value(abs_n.self_discharge_coefficient_kW_K),
            :r_state_gnft_baseline =>
                parameter_value(timeseries_to_backbone_map(result.temperatures[n])),
        ) for (n, abs_n) in result.archetype.abstract_nodes
    )
    add_relationship_parameter_values!(backbone.grid__node, gn_param_dict; merge_values=true)
    # NOTE! System link nodes handled separately via `add_system_link_node_parameters`.

    # `grid__node__boundary` based on the maximum and minimum temperatures.
    gnb_param_dict = Dict(
        (grid=g_map[n], node=n_map[n], boundary=backbone.boundary(:upwardLimit)) => Dict(
            :constant => parameter_value(abs_n.maximum_temperature_K),
            :useConstant => parameter_value(true),
        ) for (n, abs_n) in result.archetype.abstract_nodes
    )
    merge!(
        gnb_param_dict,
        Dict(
            (
                grid=g_map[n],
                node=n_map[n],
                boundary=backbone.boundary(:downwardLimit),
            ) => Dict(
                :constant => parameter_value(abs_n.minimum_temperature_K),
                :useConstant => parameter_value(true),
            ) for (n, abs_n) in result.archetype.abstract_nodes
        ),
    )
    merge!(
        gnb_param_dict,
        Dict(
            (grid=g_map[n], node=n_map[n], boundary=backbone.boundary(:reference)) =>
                Dict(
                    :constant => parameter_value(result.initial_temperatures[n]),
                    :useConstant => parameter_value(true),
                ) for (n, abs_n) in result.archetype.abstract_nodes
        ),
    )
    add_relationship_parameter_values!(backbone.grid__node__boundary, gnb_param_dict; merge_values=true)

    # `grid__node__node` based on the heat transfer coefficients
    gnn_param_dict = Dict(
        (grid=g_map[n1], node1=n_map[n1], node2=n_map[n2]) =>
            Dict(:diffCoeff => parameter_value(val)) for
        (n1, abs_n1) in result.archetype.abstract_nodes for
        (n2, val) in abs_n1.heat_transfer_coefficients_kW_K
    )
    add_relationship_parameter_values!(backbone.grid__node__node, gnn_param_dict; merge_values=true)

    # `grid__node__unit__io` based on maximum flows
    gnuio_param_dict = Dict(
        (grid=g_map[n], node=n_map[n], unit=u_map[p], io=io_map[d]) => Dict(
            :capacity => parameter_value(abs(val)),
            :conversionCoeff => parameter_value(val / abs(val)),
            :unitSize => parameter_value(abs(val) / abs_p.number_of_processes),
        ) for (p, abs_p) in result.archetype.abstract_processes for
        ((d, n), val) in abs_p.maximum_flows
    )
    add_relationship_parameter_values!(backbone.grid__node__unit__io, gnuio_param_dict; merge_values=true)

    # `unit__unittype` simply attaching `HVAC` to every unit.
    add_object!(backbone.unittype, Object(:HVAC, :unittype)) # First we need to create the new unittype.
    add_relationships!(
        backbone.unit__unittype,
        [(unit=u, unittype=backbone.unittype(:HVAC)) for u in values(u_map)],
    )
end


"""
    add_system_link_node_parameters!(
        backbone::BackboneInput,
        results::Dict{Object,ArchetypeBuildingResults};
        mod::Module = @__MODULE__,
    )

Add system link node parameters into [`BackboneInput`](@ref).
"""
function add_system_link_node_parameters!(
    backbone::BackboneInput,
    results::Dict{Object,ArchetypeBuildingResults};
    mod::Module=@__MODULE__
)
    # Initialize parameter dicts for looping
    node_balance_dict = Dict()
    influx_building_baseline_consumption = Dict()

    # Loop over the results to collect them appropriately
    for r in values(results)
        # Identify system link nodes
        system_link_nodes = mod.building_archetype__system_link_node(
            building_archetype=r.archetype.archetype,
        )
        # Create grid and node mappings for system link nodes
        g_map = Dict(
            sys_node => backbone.grid(
                mod.grid_name(
                    building_archetype=r.archetype.archetype,
                    building_node=sys_node,
                ),
            ) for sys_node in system_link_nodes
        )
        n_map = Dict(
            sys_node => backbone.node(
                mod.node_name(
                    building_archetype=r.archetype.archetype,
                    building_node=sys_node,
                ),
            ) for sys_node in system_link_nodes
        )
        # Set node balance and merge with existing values.
        merge!(
            node_balance_dict,
            Dict((grid=g_map[n], node=n_map[n]) => true for n in system_link_nodes),
        )
        # Calculate total HVAC consumption and add it to existing values.
        merge!(
            +,
            influx_building_baseline_consumption,
            Dict(
                (grid=g_map[n], node=n_map[n]) =>
                    -sum( # Backbone ts_influx is negative for consumption!
                        get(r.hvac_consumption, p, 0.0) for
                        p in mod.building_process__direction__building_node(
                            direction=mod.direction(:from_node),
                            building_node=n,
                        )
                    ) for n in system_link_nodes
            ),
        )
    end

    # Form and add the system link node parameters to the backbone inputs.
    gn_param_dict = Dict(
        k => Dict(
            :nodeBalance => parameter_value(get(node_balance_dict, k, nothing)),
            :influx_building_baseline_consumption => parameter_value(
                timeseries_to_backbone_map(
                    get(influx_building_baseline_consumption, k, nothing),
                ),
            ),
        ) for k in keys(node_balance_dict)
    )
    add_relationship_parameter_values!(backbone.grid__node, gn_param_dict; merge_values=true)
end


"""
    timeseries_to_backbone_map(ts::TimeSeries)

Convert `TimeSeries` into Backbone input data `Map` corresponding to time-varying data.
"""
function timeseries_to_backbone_map(ts::TimeSeries)
    Map(
        [:f00],
        [
            Map(
                [
                    Symbol("t" * string(i; pad=6)) for
                    (i, timestamp) in enumerate(ts.indexes)
                ],
                ts.values,
            ),
        ],
    )
end


"""
    timeseries_to_backbone_map(x::Real)

Convert `Real` into Backbone input data `Map`.
"""
function timeseries_to_backbone_map(x::Real)
    Map(
        [:f00],
        [
            Map(
                [:t000000],
                [x]
            )
        ]
    )
end