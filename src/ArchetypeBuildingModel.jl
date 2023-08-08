module ArchetypeBuildingModel

using Dates
using SpineInterface
import SpineInterface._check
using Test
using PyCall
using Interpolations

# Define the set of possible solar directions, currently hard-coded.
const solar_directions = [:horizontal, :north, :east, :south, :west]

include("types.jl")
include("base.jl")
include("process_envelope.jl")
include("process_loads.jl")
include("process_node.jl")
include("process_abstract_node.jl")
include("process_scope.jl")
include("process_system.jl")
include("process_abstract_system.jl")
include("process_weather.jl")
include("tests.jl")
include("create_backbone_input.jl")
include("create_spineopt_input.jl")
include("create_generic_input.jl")
include("solve_demand.jl")
include("solve_consumption.jl")
include("export.jl")
include("main.jl")
include("util.jl")

# Exports for the `process_archetype_buildings.jl` main program
export using_spinedb,
    run_input_data_tests,
    archetype_building_processing,
    solve_archetype_building_hvac_demand,
    initialize_result_classes!,
    add_results!,
    import_data,
    write_to_url,
    SpineOptInput,
    BackboneInput,
    GenericInput
# Exports for documentation
export BuildingNodeData, AbstractNode, AbstractProcess
# Exports for `testscript.jl`
export run_parameter_tests,
    run_object_class_tests,
    run_structure_type_tests,
    ScopeData,
    create_building_weather,
    add_object_parameter_values!,
    add_relationships!,
    WeatherData,
    EnvelopeData,
    LoadsData,
    create_building_node_network,
    create_abstract_node_network,
    BuildingProcessData,
    ArchetypeBuilding,
    ArchetypeBuildingResults

end # module
