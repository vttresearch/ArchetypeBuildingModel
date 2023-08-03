#=
    create_generic_input.jl

Create ArchetypeBuildingModel.jl structure for Spine Datastores,
mostly for debugging purposes.
=#

"""
    GenericInput

Create and store the ArchetypeBuildingModel.jl structure for Spine Data Stores.

Contains the following fields:
- `building_archetype::ObjectClass`: The original `building_archetype` object.
- `scope_data::ObjectClass`: The processed [`ScopeData`](@ref).
- `envelope_data::ObjectClass`: The processed [`EnvelopeData`](@ref).
- `building_node_data::ObjectClass`: The processed [`BuildingNodeData`](@ref).
- `building_process_data::ObjectClass`: The processed [`BuildingProcessData`](@ref).
- `loads_data::ObjectClass`: The processed [`LoadsData`](@ref).
- `abstract_node_data::ObjectClass`: The processed [`AbstractNode`](@ref).
- `abstract_process_data::ObjectClass`: The processed [`AbstractProcess`](@ref).
- `building_archetype__scope_data::RelationshipClass`: Link archetype to its [`ScopeData`](@ref).
- `building_archetype__envelope_data::RelationshipClass`: Link archetype to its [`EnvelopeData`](@ref).
- `building_archetype__building_node_data::RelationshipClass`: Link archetype to its [`BuildingNodeData`](@ref).
- `building_archetype__building_process_data::RelationshipClass`: Link archetype to its [`BuildingProcessData`](@ref).
- `building_archetype__loads_data::RelationshipClass`: Link archetype to its [`LoadsData`](@ref).
- `building_archetype__abstract_node_data::RelationshipClass`: Link archetype to its [`AbstractNode`](@ref).
- `building_archetype__abstract_process_data::RelationshipClass`: Link archetype to its [`AbstractProcess`](@ref).
"""
struct GenericInput <: ModelInput
    building_archetype::ObjectClass
    scope_data::ObjectClass
    envelope_data::ObjectClass
    building_node_data::ObjectClass
    building_process_data::ObjectClass
    loads_data::ObjectClass
    abstract_node_data::ObjectClass
    abstract_process_data::ObjectClass
    building_archetype__scope_data::RelationshipClass
    building_archetype__envelope_data::RelationshipClass
    building_archetype__building_node_data::RelationshipClass
    building_archetype__building_process_data::RelationshipClass
    building_archetype__loads_data::RelationshipClass
    building_archetype__abstract_node_data::RelationshipClass
    building_archetype__abstract_process_data::RelationshipClass
    function GenericInput()
        building_archetype = ObjectClass(:building_archetype, Array{ObjectLike,1}())
        scope_data = ObjectClass(:scope_data, Array{ObjectLike,1}())
        envelope_data = ObjectClass(:envelope_data, Array{ObjectLike,1}())
        building_node_data = ObjectClass(:building_node_data, Array{ObjectLike,1}())
        building_process_data = ObjectClass(:building_process_data, Array{ObjectLike,1}())
        loads_data = ObjectClass(:loads_data, Array{ObjectLike,1}())
        abstract_node_data = ObjectClass(:abstract_node_data, Array{ObjectLike,1}())
        abstract_process_data = ObjectClass(:abstract_process_data, Array{ObjectLike,1}())
        building_archetype__scope_data = RelationshipClass(
            :building_archetype__scope_data,
            [:building_archetype, :scope_data],
            Array{RelationshipLike,1}(),
        )
        building_archetype__envelope_data = RelationshipClass(
            :building_archetype__envelope_data,
            [:building_archetype, :envelope_data],
            Array{RelationshipLike,1}()
        )
        building_archetype__building_node_data = RelationshipClass(
            :building_archetype__building_node_data,
            [:building_archetype, :building_node_data],
            Array{RelationshipLike,1}()
        )
        building_archetype__building_process_data::RelationshipClass
        building_archetype__loads_data::RelationshipClass
        building_archetype__abstract_node_data::RelationshipClass
        building_archetype__abstract_process_data::RelationshipClass
    end
end