#=
    process_archetype_buildings.jl

The main program for running the archetype building model via Spine Toolbox.
=#

using Pkg
Pkg.activate(@__DIR__)
using ArchetypeBuildingModel
using Test

# Check that the necessary input arguments are provided
if length(ARGS) < 1
    @error """
    `process_archetype_buildings` requires at least the following input arguments:
    1. The url to a Spine Datastore containing the required input data and archetype building definitions.

    Furthermore, the following optional keyword arguments can be provided:
    2. `-spineopt <url>`, the url to a Spine Datastore where the produced SpineOpt input data should be written.
    3. `-backbone <url>`, the url to a Spine Datastore where the produced Backbone input data should be written.
    4. `-import_weather <false>`, controls whether auto-generated `building_weather` are imported into the DB.
    5. `-save_layouts <false>`, controls whether auto-generated `building_weather` layouts are saved as images.
    """
else
    # Process command line arguments
    url_in = popfirst!(ARGS)
    kws = Dict(ARGS[i] => get(ARGS, i + 1, nothing) for i = 1:2:length(ARGS))
    spineopt_url = get(kws, "-spineopt", nothing)
    backbone_url = get(kws, "-backbone", nothing)
    import_weather = lowercase(get(kws, "-import_weather", "false")) == "true"
    save_layouts = lowercase(get(kws, "-save_layouts", "false")) == "true"

    # Open input database and run tests.
    @info "Opening input datastore at `$(url_in)`..."
    @time using_spinedb(url_in, Main)

    # Run input data tests
    run_input_data_tests()

    # Process ScopeData and WeatherData, and create the ArchetypeBuildings
    scope_data_dictionary, weather_data_dictionary, archetype_dictionary =
        archetype_building_processing(url_in, import_weather)

    # Heating/cooling demand calculations.
    archetype_results_dictionary =
        solve_archetype_building_hvac_demand(archetype_dictionary; free_dynamics = false)

    # Write the results back into the input datastore
    results__building_archetype__building_node,
    results__building_archetype__building_process = initialize_result_relationship_classes()
    add_results!(
        results__building_archetype__building_node,
        results__building_archetype__building_process,
        archetype_results_dictionary,
    )
    @info "Importing `ArchetypeBuildingResults`..."
    @time import_data(
        url_in,
        [
            results__building_archetype__building_node,
            results__building_archetype__building_process,
        ],
        "Importing `ArchetypeBuildingResults`.",
    )

    # Process SpineOpt input data if requested.
    if !isnothing(spineopt_url)
        @info "Processing and writing SpineOpt input data into `$(spineopt_url)`..."
        @time write_to_url(
            String(spineopt_url),
            SpineOptInput(archetype_dictionary, archetype_results_dictionary),
        )
    end

    # Process Backbone input data if requested.
    if !isnothing(backbone_url)
        @info "Processing and writing Backbone input data into `$(backbone_url)`..."
        @time write_to_url(
            String(backbone_url),
            BackboneInput(archetype_dictionary, archetype_results_dictionary),
        )
    end

    @info """
    All done!
    You can access the `ArchetypeBuilding` data in the `archetype_dictionary`,
    and the `ArchetypeBuildingResults` in the `archetype_results_dictionary`.
    `ScopeData` and `WeatherData` are also available in the `scope_data_dictionary`
    and `weather_data_dictionary` respectively.
    """
end
