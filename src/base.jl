#=
    base.jl

Contains functions affecting Julia `Base` module.
=#


# Define custom `show` for `ArchetypeBuilding` to avoid overloading REPL.
function Base.show(io::IO, archetype::ArchetypeBuilding)
    println(
        io,
        """
        ArchetypeBuilding
            archetype:  $(archetype.archetype)
            scope:      $(archetype.scope)
            fabrics:    $(archetype.fabrics)
            systems:    $(archetype.systems)
            loads:      $(archetype.loads)
            weather:    $(archetype.weather)
        """,
    )
end


# Define custom `show` for `BuildingDataType`
function Base.show(io::IO, building_data::BuildingDataType)
    str = ["$(typeof(building_data))\n"]
    for field in fieldnames(typeof(building_data))
        push!(str, "    $(field):\t$(getfield(building_data, field))\n")
    end
    println(io, join(str))
end
