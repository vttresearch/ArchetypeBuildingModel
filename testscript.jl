#=
    testscript.jl

A makeshift script for testing the functionality of `ArchetypeBuildingModel.jl`.

Contains a lot of old code snippets that might not work without changes to the
module. Included primarily to give some idea how the module operates under the
hood.
=#

using Pkg
Pkg.activate("test")

using Revise
using Test
using Plots
using ArchetypeBuildingModel
m = Module()

# Open database

# Provide the url for a datastore containing the required raw input data and the archetype building definitions.
url = "sqlite:///C:\\_SPINEPROJECTS\\mopo_ambience_test\\.spinetoolbox\\items\\abm_data\\ABM_data.sqlite"

# Output url
#output_url = <ADD OUTPUT URL IF DESIRED>

@info "Opening database..."
@time using_spinedb(url, m)
realization = :realization # Realization data scenario name from stochastic input


## Run tests

lmt = Inf
@info "Running input data and definition tests..."
@testset begin
    @time run_parameter_tests(m; limit=lmt)
    @time run_object_class_tests(m; limit=lmt)
    @time run_structure_type_tests(m)
end

## Test creating the `ScopeData` types automatically
#=
@info "Processing the `ScopeData` objects for the test `building_scope` objects..."
@time scope_data_final = Dict(scope => ScopeData(scope) for scope in building_scope())
=#


## Test autocreation of `WeatherData`
#=
@info "Testing automatic creation of `building_weather`..."
@time auto_weather = Dict(
    archetype => create_building_weather(
        archetype,
        scope_data_final[first(
            building_archetype__building_scope(building_archetype = archetype),
        )],
    ) for archetype in building_archetype()
)
# Add the created `building_weather` objects into the in-memory DB.
for (archetype, (bw, bw_params)) in auto_weather
    add_object_parameter_values!(building_weather, bw_params)
    add_relationships!(
        building_archetype__building_weather,
        [(building_archetype = archetype, building_weather = bw)],
    )
end
=#


## Test creating `WeatherData`
#=
@info "Processing `WeatherData` objects..."
@time weather_data = Dict(weather => WeatherData(weather) for weather in building_weather())

temperature_plot = plot(first(weather_data)[2].ambient_temperature_K.values)
temperature_plot = plot!(first(weather_data)[2].ground_temperature_K.values)
display(temperature_plot)
=#


## Test creating `EnvelopeData`
#=
@info "Processing the `EnvelopeData` objects for the test `building_archetype` objects..."
@time envelope_data = Dict(
    archetype => EnvelopeData(
        archetype,
        scope_data_final[first(
            building_archetype__building_scope(building_archetype = archetype),
        )],
    ) for archetype in building_archetype()
)
=#


## Test creating `LoadsData`
#=
@info "Processing the `LoadsData` objects for the test `building_archetype` objects..."
@time loads_data = Dict(
    archetype => LoadsData(
        archetype,
        scope_data_final[first(
            building_archetype__building_scope(building_archetype = archetype),
        )],
        envelope_data[archetype],
        weather_data[first(
            building_archetype__building_weather(building_archetype = archetype),
        )],
    ) for archetype in building_archetype()
)
=#


## Test creating `BuildingNodeNetwork` and `BuildingNodeData`
#=
@info "Processing `BuildingNodeNetwork` and `BuildingNodeData` for the test `building_archetype` objects..."
@time building_node_network = Dict(
    archetype => create_building_node_network(
        archetype,
        first(building_archetype__building_fabrics(building_archetype = archetype)),
        first(building_archetype__building_systems(building_archetype = archetype)),
        scope_data_final[first(
            building_archetype__building_scope(building_archetype = archetype),
        )],
        envelope_data[archetype],
        loads_data[archetype],
    ) for archetype in building_archetype()
)
=#


## Test creating `AbstractNodeNetwork` and `AbstractNode`
#=
@info "Processing `AbstractNodeNetwork` and `AbstractNode`..."
@time abstract_node_network = Dict(
    archetype => create_abstract_node_network(
        building_node_network[archetype],
        weather_data[first(
            building_archetype__building_weather(building_archetype = archetype),
        )],
    ) for archetype in building_archetype()
)
=#


## Test creating `BuildingProcessData`
#=
@info "Processing the `BuildingProcessData` objects for the test `building_archetype` objects..."
@time building_process_data = Dict(
    (archetype, process) => BuildingProcessData(
        archetype,
        process,
        scope_data_final[first(
            building_archetype__building_scope(building_archetype = archetype),
        )],
        weather_data[first(
            building_archetype__building_weather(building_archetype = archetype),
        )],
    ) for archetype in building_archetype() for
    process in building_systems__building_process(
        building_systems = first(
            building_archetype__building_systems(building_archetype = archetype),
        ),
    )
)
=#


## Test time-dependent COPs, ONLY WORKS FOR VERY SPECIFIC INPUT DATA!
#=
archetype = first(building_archetype())
scopedata = first(scope_data_final)[2]
weatherdata = first(weather_data)[2]

# Test different heat pumps.
a2ahp = BuildingProcessData(archetype, building_process(:A2AHP), scopedata, weatherdata)
g2whp = BuildingProcessData(archetype, building_process(:G2WHP), scopedata, weatherdata)

plot(
    weatherdata.ambient_temperature_K.indexes,
    weatherdata.ambient_temperature_K.values .- 273.15,
)
plot!(a2ahp.coefficient_of_performance.indexes, a2ahp.coefficient_of_performance.values)
plot!(g2whp.coefficient_of_performance.indexes, g2whp.coefficient_of_performance.values)
=#


## Test creating `ArchetypeBuilding`s

@info "Processing `ArchetypeBuilding` objects..."
@time archetype_dictionary = Dict(
    archetype => ArchetypeBuilding(archetype; mod=m, realization=realization) for
    archetype in m.building_archetype()
)


## Test heating/cooling demand calculations.

@info "Creating `ArchetypeBuildingResults`..."
@time archetype_results = Dict(
    archetype => ArchetypeBuildingResults(
        val;
        free_dynamics=false,
        mod=m,
        realization=realization
    ) for (archetype, val) in archetype_dictionary
)


## Test creating the results database structures.

@info "Initializing result classes..."
@time results__building_archetype__building_node,
results__building_archetype__building_process,
results__system_link_node = initialize_result_classes!(m)
@info "Adding results..."
@time add_results!(
    results__building_archetype__building_node,
    results__building_archetype__building_process,
    results__system_link_node,
    archetype_results;
    mod=m
)


## Test creating and writing SpineOpt input

@info "Creating `SpineOptInput`..."
@time spineopt = SpineOptInput(archetype_results; mod=m)
#write_to_url(output_url, spineopt)


## Test creating and writing Backbone input

@info "Creating `BackboneInput`..."
@time backbone = BackboneInput(archetype_results; mod=m)
#write_to_url(output_url, backbone)


## Plot diagnostics.

results = first(values(archetype_results))

weather_plt = plot(; title="Ambient temperatures in [C]")
plot!(
    weather_plt,
    keys(results.archetype.weather_data.ambient_temperature_K),
    values(results.archetype.weather_data.ambient_temperature_K) .- 273.15;
    label="Ambient"
)
plot!(
    weather_plt,
    keys(results.archetype.weather_data.ground_temperature_K),
    values(results.archetype.weather_data.ground_temperature_K) .- 273.15;
    label="Ground"
)
display(weather_plt)

temp_plt = plot(; title="Node temperatures in [C]")
for (n, ts) in results.temperatures
    plot!(temp_plt, keys(ts), values(ts) .- 273.15, label=string(n))
end
display(temp_plt)

hvac_plt = plot(; title="Heating/cooling demand in [W]")
for (n, ts) in results.hvac_demand
    plot!(hvac_plt, keys(ts), values(ts), label=string(n))
end
display(hvac_plt)

process_plt = plot(; title="HVAC consumption in [MW]")
for (p, ts) in results.hvac_consumption
    plot!(process_plt, keys(ts), values(ts), label=string(p))
end
display(process_plt)
