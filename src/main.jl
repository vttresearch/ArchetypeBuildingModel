#=
    main.jl

Contains functions for the main program `process_archetype_buildings.jl`.
=#

"""
    run_input_data_tests(mod::Module = @__MODULE__)

Runs input data tests for the Datastore loaded to module `mod`, `@__MODULE__` by default.

Essentially performs the following steps:
1. Call [`run_object_class_tests`](@ref)
2. Call [`run_parameter_tests`](@ref)
3. Call [`run_structure_type_tests`](@ref)
"""
function run_input_data_tests(mod::Module = @__MODULE__)
    @time @testset "Datastore tests" begin
        run_object_class_tests(mod)
        run_parameter_tests(mod)
        run_structure_type_tests(mod)
    end
end


"""
    archetype_building_processing(
        weather_url::String,
        save_layouts::Bool;
        weather_data_dictionary::Union{Nothing,Dict{Object,WeatherData}} = nothing,
        mod::Module = @__MODULE__,
        realization::Symbol = :realization,
    )

Process the [`ScopeData`](@ref), [`WeatherData`](@ref), and [`ArchetypeBuilding`](@ref) objects.

Essentially, processes all the necessary information for [`ArchetypeBuilding`](@ref)
creation, and returns the `scope_data_dictionary`, `weather_data_dictionary`,
and `archetype_dictionary` for examining the processed data.
Any automatically generated [building_weather](@ref)
objects will be imported back into the database at `weather_url`.
If `save_layouts == true`, diagnostic figures of the layouts are saved into `figs/`.
The `weather_data_dictionary` keyword can be used to bypass weather data processing
if a pre-existing dictionary is provided.
The `mod` keyword changes from which Module data is accessed from, `@__MODULE__` by default.
The `realization` keyword is used to indicate the true data from potentially
stochastic input.

This function performs the following steps:
1. Construct the [`ScopeData`](@ref) for each defined [building\\_archetype\\_\\_building_scope](@ref), and store in the `scope_data_dictionary`.
2. Try to construct the [`WeatherData`](@ref) for each defined [building\\_archetype\\_\\_building_weather](@ref), and attempt automatic weather processing using [ArchetypeBuildingWeather.py](@ref) if no definition found. Results stored in `weather_data_dictionary`.
3. Use the `scope_data_dictionary` and `weather_data_dictionary` to construct the [`ArchetypeBuilding`](@ref) for all defined archetypes, and store them in `archetype_dictionary`.
4. Return `scope_data_dictionary`, `weather_data_dictionary`, and `archetype_dictionary`.
"""
function archetype_building_processing(
    weather_url::String,
    save_layouts::Bool;
    weather_data_dictionary::Union{Nothing,Dict{Object,WeatherData}} = nothing,
    mod::Module = @__MODULE__,
    realization::Symbol = :realization,
)
    # Process relevant `ScopeData` objects.
    @info "Processing `building_scope` objects into `ScopeData` for `scope_data_dictionary`..."
    @time scope_data_dictionary = Dict(
        archetype => ScopeData(scope; mod = mod) for
        (archetype, scope) in mod.building_archetype__building_scope()
    )

    # Process relevant `WeatherData` objects.
    if isnothing(weather_data_dictionary)
        archetypes_missing_weather = setdiff(
            mod.building_archetype(),
            getfield.(mod.building_archetype__building_weather(), :building_archetype),
        )
        if !isempty(archetypes_missing_weather)
            @warn """
            Creating missing `building_weather` objects using `ArchetypeBuildingWeather.py`.
            CAUTION! This might take a while if the weather data isn't readily available.
            """
            @time for archetype in archetypes_missing_weather
                bw, bw_params = create_building_weather(
                    archetype,
                    scope_data_dictionary[archetype];
                    save_layouts = save_layouts,
                )
                add_object_parameter_values!(mod.building_weather, bw_params)
                add_relationships!(
                    mod.building_archetype__building_weather,
                    [(building_archetype = archetype, building_weather = bw)],
                )
            end
            @info "Importing auto-generated `building_weather` into `$(weather_url)`..."
            @time import_data(
                weather_url,
                [mod.building_weather, mod.building_archetype__building_weather],
                "Autogenerated `building_weather` and `building_archetype__building_weather`",
            )
        end
        @info "Processing `building_weather` objects into `WeatherData` for `weather_data_dictionary`..."
        @time weather_data_dictionary = Dict(
            archetype => WeatherData(weather; mod = mod, realization = realization) for
            (archetype, weather) in mod.building_archetype__building_weather()
        )
    else
        @info "Using given `weather_data_dictionary`."
    end

    # Process `ArchetypeBuilding` objects.
    @info "Processing `building_archetype` objects into `ArchetypeBuilding` for `archetype_dictionary`..."
    @time archetype_dictionary = Dict(
        archetype => ArchetypeBuilding(
            archetype,
            scope_data_dictionary[archetype],
            weather_data_dictionary[archetype];
            mod = mod,
        ) for archetype in mod.building_archetype()
    )

    # Return the dictionaries of interest
    return scope_data_dictionary, weather_data_dictionary, archetype_dictionary
end


"""
    solve_archetype_building_hvac_demand(
        archetype_dictionary::Dict{Object,ArchetypeBuilding};
        free_dynamics::Bool = false,
        initial_temperatures::Dict{Object,Dict{Object,Float64}} = Dict{
            Object,
            Dict{Object,Float64}
        }(),
        realization::Symbol = :realization,
    )

Solve the [`ArchetypeBuilding`](@ref) heating and cooling demand.

The `free_dynamics` keyword can be used to ignore node temperature limits,
while the `initial_temperatures` keyword can be used to set desired initial
temperatures for the nodes. The `realization` keyword is used to denote
the true data from potentially stochastic input.

Essentially, performs the following steps:
1. Create the `archetype_results_dictionary` by constructing the [`ArchetypeBuildingResults`](@ref) for each entry in the `archetype_dictionary`.
2. Create the `results__building_archetype__building_node` `RelationshipClass` for storing temperature results.
3. Create the `results__building_archetype__building_process` `RelationshipClass` for storing HVAC results.
4. Return the `archetype_results_dictionary`, as well as the created `RelationshipClasses`.
"""
function solve_archetype_building_hvac_demand(
    archetype_dictionary::Dict{Object,ArchetypeBuilding};
    free_dynamics::Bool = false,
    initial_temperatures::Dict{Object,Dict{Object,Float64}} = Dict{
        Object,
        Dict{Object,Float64},
    }(),
    mod::Module = @__MODULE__,
    realization::Symbol = :realization,
)
    # Heating/cooling demand calculations.
    @info "Calculating heating/cooling demand..."
    @time archetype_results_dictionary = Dict(
        archetype => ArchetypeBuildingResults(
            archetype_building;
            free_dynamics = free_dynamics,
            initial_temperatures = get(initial_temperatures, archetype, nothing),
            mod = mod,
            realization = realization,
        ) for (archetype, archetype_building) in archetype_dictionary
    )

    # Return the results dictionary
    return archetype_results_dictionary
end


"""
    initialize_result_classes!(mod::Module)

Initialize `RelationshipClass`es and `ObjectClass`es for storing heating and HVAC demand results in `mod`.

Note that this function modifies `mod` directly!
"""
function initialize_result_classes!(mod::Module)
    # Initialize archetype node results
    results__building_archetype__building_node = RelationshipClass(
        :results__building_archetype__building_node,
        [:building_archetype, :building_node],
        Array{RelationshipLike,1}(),
        Dict(),
        Dict(
            param => parameter_value(nothing) for
            param in [:initial_temperature_K, :temperature_K, :hvac_demand_W]
        ),
    )
    # Create the associated parameters
    initial_temperature_K =
        Parameter(:initial_temperature_K, [results__building_archetype__building_node])
    temperature_K = Parameter(:temperature_K, [results__building_archetype__building_node])
    hvac_demand_W = Parameter(:hvac_demand_W, [results__building_archetype__building_node])

    # Initialize process results
    results__building_archetype__building_process = RelationshipClass(
        :results__building_archetype__building_process,
        [:building_archetype, :building_process],
        Array{RelationshipLike,1}(),
        Dict(),
        Dict(:hvac_consumption_MW => parameter_value(nothing)),
    )
    # Create the associated parameter
    hvac_consumption_MW =
        Parameter(:hvac_consumption_MW, [results__building_archetype__building_process])

    # Initialize system link node results
    results__system_link_node = ObjectClass(
        :results__system_link_node,
        Array{ObjectLike,1}(),
        Dict(),
        Dict(:total_consumption_MW => parameter_value(nothing)),
    )
    # Create the assosicated parameter
    total_consumption_MW = Parameter(:total_consumption_MW, [results__system_link_node])

    # Evaluate the relationship classes and parameters to the desired module.
    @eval mod begin
        results__building_archetype__building_node =
            $results__building_archetype__building_node
        results__building_archetype__building_process =
            $results__building_archetype__building_process
        results__system_link_node = $results__system_link_node
        initial_temperature_K = $initial_temperature_K
        temperature_K = $temperature_K
        hvac_demand_W = $hvac_demand_W
        hvac_consumption_MW = $hvac_consumption_MW
        total_consumption_MW = $total_consumption_MW
    end

    # Return the handles for the relationship classes for future reference.
    return results__building_archetype__building_node,
    results__building_archetype__building_process,
    results__system_link_node
end


"""
    add_results!(
        results__building_archetype__building_node::RelationshipClass,
        results__building_archetype__building_process::RelationshipClass,
        results__system_link_node::ObjectClass,
        results_dictionary::Dict{Object,ArchetypeBuildingResults};
        mod::Module = @__MODULE__
    )

    Add results from `results_dictionary` into the result `RelationshipClass`es.

    NOTE! The `mod` keyword changes from which Module data is accessed from,
    `@__MODULE__` by default.
"""
function add_results!(
    results__building_archetype__building_node::RelationshipClass,
    results__building_archetype__building_process::RelationshipClass,
    results__system_link_node::ObjectClass,
    results_dictionary::Dict{Object,ArchetypeBuildingResults};
    mod::Module = @__MODULE__,
)
    # Collect `ArchetypeBuildingResults`
    results = values(results_dictionary)

    # Add `results__building_archetype__building_node` results.
    add_relationship_parameter_values!(
        results__building_archetype__building_node,
        Dict(
            (building_archetype = r.archetype.archetype, building_node = node) => Dict(
                :initial_temperature_K => parameter_value(r.initial_temperatures[node]),
                :temperature_K => parameter_value(r.temperatures[node]),
                :hvac_demand_W => parameter_value(r.hvac_demand[node]),
            ) for r in results for node in keys(r.temperatures)
        ),
    )

    # Add `results__building_archetype__building_process` results.
    add_relationship_parameter_values!(
        results__building_archetype__building_process,
        Dict(
            (building_archetype = r.archetype.archetype, building_process = process) =>
                Dict(:hvac_consumption_MW => parameter_value(r.hvac_consumption[process]))
            for r in results for process in keys(r.hvac_consumption)
        ),
    )

    # Add `results__system_link_node` results.
    total_cons_MW = merge(+, getfield.(results, :hvac_consumption)...)
    add_object_parameter_values!(
        results__system_link_node,
        Dict(
            sys_link_n => Dict(
                :total_consumption_MW => parameter_value(
                    sum(
                        get(total_cons_MW, p, 0.0) for
                        p in mod.building_process__direction__building_node(
                            direction = mod.direction(:from_node),
                            building_node = sys_link_n,
                        )
                    ),
                ),
            ) for sys_link_n in mod.building_archetype__system_link_node(
                building_archetype = mod.building_archetype(),
            )
        ),
    )

    # Return the results of interest
    return results__building_archetype__building_node,
    results__building_archetype__building_process,
    results__system_link_node
end