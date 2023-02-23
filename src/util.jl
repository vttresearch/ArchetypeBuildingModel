#=
    util.jl

Contains various utility functions for convenience.
=#


"""
    collect_leaf_values(x::SpineDataType)

Returns an `Array` of the leaf values contained in any `SpineDataType`.
"""
function collect_leaf_values(x::SpineDataType)
    return collect(values(x))
end
function collect_leaf_values(x::Map)
    return vcat(collect_leaf_values.(values(x))...)
end
