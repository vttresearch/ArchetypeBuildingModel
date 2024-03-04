#=
    testscript.jl

A makeshift script for testing the functionality of `ArBuMo.jl`.

Contains a lot of old code snippets that might not work without changes to the
module. Included primarily to give some idea how the module operates under the
hood.
=#

using Pkg
Pkg.activate("test")

using Test
using Plots
using ArBuMo
m = Module()

# Open database

# Provide the url for a datastore containing the required raw input data and the archetype building definitions.
url = "sqlite:///C:\\_SPINEPROJECTS\\SpineOpt_PED_demo_fluid\\.spinetoolbox\\data_and_definitions.sqlite"

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

@info "Processing the `ScopeData` objects for the test `building_scope` objects..."
@time scope_data_final = Dict(scope => ScopeData(scope; mod=m) for scope in m.building_scope())


## Test creating `EnvelopeData`

@info "Processing the `EnvelopeData` objects for the test `building_archetype` objects..."
@time envelope_data = Dict(
    archetype => EnvelopeData(
        archetype,
        scope_data_final[only(
            m.building_archetype__building_scope(building_archetype=archetype),
        )];
        mod=m,
    ) for archetype in m.building_archetype()
)


## Test creating `LoadsData`

@info "Processing the `LoadsData` objects for the test `building_archetype` objects..."
@time loads_data = Dict(
    archetype => LoadsData(
        archetype,
        scope_data_final[only(
            m.building_archetype__building_scope(building_archetype=archetype),
        )],
        envelope_data[archetype];
        mod=m,
    ) for archetype in m.building_archetype()
)


## Test creating `BuildingNodeNetwork` and `BuildingNodeData`

@info "Processing `BuildingNodeNetwork` and `BuildingNodeData` for the test `building_archetype` objects..."
@time building_node_network = Dict(
    archetype => create_building_node_network(
        archetype,
        only(m.building_archetype__building_fabrics(building_archetype=archetype)),
        only(m.building_archetype__building_systems(building_archetype=archetype)),
        scope_data_final[only(
            m.building_archetype__building_scope(building_archetype=archetype),
        )],
        envelope_data[archetype],
        loads_data[archetype];
        mod=m
    ) for archetype in m.building_archetype()
)


## Test autocreation of `WeatherData` 

@info "Testing automatic creation of `building_weather`..."
@time weather_data = Dict(
    archetype => WeatherData(
        archetype,
        scope_data_final[only(
            m.building_archetype__building_scope(building_archetype=archetype),
        )],
        envelope_data[archetype],
        building_node_network[archetype],
        loads_data[archetype];
        mod=m,
    ) for archetype in m.building_archetype()
)


## Test creating `BuildingProcessData`

@info "Processing the `BuildingProcessData` objects for the test `building_archetype` objects..."
@time building_process_data = Dict(
    (archetype, process) => BuildingProcessData(
        archetype,
        process,
        scope_data_final[only(
            m.building_archetype__building_scope(building_archetype=archetype),
        )],
        weather_data[archetype];
        mod=m,
    ) for archetype in m.building_archetype() for
    process in m.building_systems__building_process(
        building_systems=only(
            m.building_archetype__building_systems(building_archetype=archetype),
        ),
    )
)


## Test creating `AbstractNodeNetwork` and `AbstractNode`

@info "Processing `AbstractNodeNetwork` and `AbstractNode`..."
@time abstract_node_network = Dict(
    archetype => create_abstract_node_network(
        building_node_network[archetype]
    ) for archetype in m.building_archetype()
)


## Test creating `AbstractProcess`es.

@info "Processing `AbstractProcess`..."
@time abstract_process_data = Dict(
    (archetype, process) => AbstractProcess(process_data; mod=m)
    for ((archetype, process), process_data) in building_process_data
)


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


## Test creating generic input

@info "Creating `GenericInput`..."
@time generic = GenericInput(archetype_results; mod=m)
#@time write_to_url(output_url, generic)


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
