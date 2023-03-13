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
