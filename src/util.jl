#=
    util.jl

Contains various utility functions for convenience.
=#


"""
    collect_leaf_values(x::T) where T ∈ {SpineDataType,AbstractArray,AbstractDict}

Returns an `Array` of the leaf values contained in any `SpineDataType`.
"""
function collect_leaf_values(x::Union{Real,TimeSeries,TimePattern})
    return collect(values(x))
end
function collect_leaf_values(x::Union{Map,AbstractDict,AbstractArray})
    return vcat(collect_leaf_values.(values(x))...)
end


"""
    load_definitions_template()

Read the `archetype_definitions.json` into `Dict``.
"""
function load_definitions_template()
    JSON.parsefile(@__DIR__ * "\\..\\archetype_definitions.json")
end