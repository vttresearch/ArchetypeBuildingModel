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
- `building_nodes::ObjectClass`: The processed [`BuildingNodeData`](@ref).
- `building_processes::ObjectClass`: The processed [`BuildingProcessData`](@ref).
- `loads_data::ObjectClass`: The processed [`LoadsData`](@ref).
- `abstract_nodes::ObjectClass`: The processed [`AbstractNode`](@ref).
- `abstract_processes::ObjectClass`: The processed [`AbstractProcess`](@ref).
- `building_archetype__scope_data::RelationshipClass`: Link archetype to its [`ScopeData`](@ref).
- `building_archetype__envelope_data::RelationshipClass`: Link archetype to its [`EnvelopeData`](@ref).
- `building_archetype__building_nodes::RelationshipClass`: Link archetype to its [`BuildingNodeData`](@ref).
- `building_archetype__building_processes::RelationshipClass`: Link archetype to its [`BuildingProcessData`](@ref).
- `building_archetype__loads_data::RelationshipClass`: Link archetype to its [`LoadsData`](@ref).
- `building_archetype__abstract_nodes::RelationshipClass`: Link archetype to its [`AbstractNode`](@ref).
- `building_archetype__abstract_processes::RelationshipClass`: Link archetype to its [`AbstractProcess`](@ref).
"""
struct GenericInput <: ModelInput
    building_archetype::ObjectClass
    scope_data::ObjectClass
    envelope_data::ObjectClass
    building_nodes::ObjectClass
    building_processes::ObjectClass
    loads_data::ObjectClass
    abstract_nodes::ObjectClass
    abstract_processes::ObjectClass
    building_archetype__scope_data::RelationshipClass
    building_archetype__envelope_data::RelationshipClass
    building_archetype__building_nodes::RelationshipClass
    building_archetype__building_processes::RelationshipClass
    building_archetype__loads_data::RelationshipClass
    building_archetype__abstract_nodes::RelationshipClass
    building_archetype__abstract_processes::RelationshipClass
    function GenericInput()
        building_archetype = ObjectClass(:building_archetype, Array{ObjectLike,1}())
        scope_data = ObjectClass(:scope_data, Array{ObjectLike,1}())
        envelope_data = ObjectClass(:envelope_data, Array{ObjectLike,1}())
        building_nodes = ObjectClass(:building_nodes, Array{ObjectLike,1}())
        building_processes = ObjectClass(:building_processes, Array{ObjectLike,1}())
        loads_data = ObjectClass(:loads_data, Array{ObjectLike,1}())
        abstract_nodes = ObjectClass(:abstract_nodes, Array{ObjectLike,1}())
        abstract_processes = ObjectClass(:abstract_processes, Array{ObjectLike,1}())
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
        building_archetype__building_nodes = RelationshipClass(
            :building_archetype__building_nodes,
            [:building_archetype, :building_nodes],
            Array{RelationshipLike,1}()
        )
        building_archetype__building_processes = RelationshipClass(
            :building_archetype__building_processes,
            [:building_archetype, :building_processes],
            Array{RelationshipLike,1}()
        )
        building_archetype__loads_data = RelationshipClass(
            :building_archetype__loads_data,
            [:building_archetype, :loads_data],
            Array{RelationshipLike,1}()
        )
        building_archetype__abstract_nodes = RelationshipClass(
            :building_archetype__abstract_nodes,
            [:building_archetype, :abstract_nodes],
            Array{RelationshipLike,1}()
        )
        building_archetype__abstract_processes = RelationshipClass(
            :building_archetype__abstract_processes,
            [:building_archetype, :abstract_processes],
            Array{RelationshipLike,1}()
        )
        new(
            building_archetype,
            scope_data,
            envelope_data,
            building_nodes,
            building_processes,
            loads_data,
            abstract_nodes,
            abstract_processes,
            building_archetype__scope_data,
            building_archetype__envelope_data,
            building_archetype__building_nodes,
            building_archetype__building_processes,
            building_archetype__loads_data,
            building_archetype__abstract_nodes,
            building_archetype__abstract_processes,
        )
    end
end


"""
    GenericInput(
        archetypes::Dict{Object,ArchetypeBuilding}
    )

Create [`GenericInput`](@ref) based on a given archetype building dictionary.

Essentially, performs the following steps:
1. Initialize an empty [`GenericInput`](@ref).
2. Loop over the given `archetypes`, and [`add_archetype_to_input!`](@ref) one by one.
"""
function GenericInput(
    archetypes::Dict{Object,ArchetypeBuilding}
)
    generic = GenericInput()
    for archetype in values(archetypes)
        add_archetype_to_input!(generic, archetype)
    end
    return generic
end


"""
    add_archetype_to_input!(
        generic::GenericInput,
        archetype::ArchetypeBuilding
    )

Add [`ArchetypeBuilding`](@ref) to [`GenericInput`](@ref).

Essentially loops over the 
"""
function add_archetype_to_input!(
    generic::GenericInput,
    archetype::ArchetypeBuilding
)
    # Define the fields of ArchetypeBuilding to loop over first.
    fields = [ # These don't contain deeper levels.
        :scope_data,
        :envelope_data,
        :loads_data,
        :building_nodes,
        :building_processes,
        :abstract_nodes,
        :abstract_processes,
    ]

    # Generate the required archetype and field objects with mapping.
    ao_map = Dict(
        (archetype, field, data) => Object(
            Symbol(string(archetype.archetype, name) * "__" * string((getfield(data, 1))))
        )
        for field in fields
        for data in values(getfield(archetype, field))
    )

    # Loop over the mappings to populate GenericInput
    for ((archetype, field, data), obj) in ao_map
        # Process objects and their properties into GenericInput
        properties = fieldnames(typeof(data))[2:end] # Skip first property as it contains the object.
        param_dict = Dict(
            obj => Dict(
                key => parameter_value(getfield(data, key))
                for key in properties
            )
        )
        add_object_parameter_values!(getfield(generic, field), param_dict)
        merge!(
            getfield(getfield(generic, field), :parameter_defaults),
            Dict(
                key => parameter_value(nothing)
                for key in properties
            )
        )
        # Process relationships
        relclass = Symbol("building_archetype__" * string(field))
        add_relationships!(
            getfield(generic, relclass),
            [NamedTuple{(:building_archetype, field)}(archetype, obj)]
        )
    end

    return nothing
end


"""
    SpineInterface.parameter_value(d::Dict)

Extend `SpineInterface.parameter_value` to dictionaries.
"""
function SpineInterface.parameter_value(d::Dict)
    return parameter_value(
        Map(
            keys(d),
            parameter_value.(values(d))
        )
    )
end


"""
    SpineInterface.parameter_value(sd::StructureData)

Extend `SpineInterface.parameter_value` to [`StructureData`](@ref).
"""
function SpineInterface.parameter_value(sd::StructureData)
    keys = fieldnames(sd)[2:end] # Skip first value as it's the structure type object.
    return parameter_value(
        Map(
            keys,
            [getfield(sd, k) for k in keys]
        )
    )
end