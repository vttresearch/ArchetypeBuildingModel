#=
    solve_consumption.jl

Functions to solve the energy consumption of the `abstract_processes`,
not just the HVAC demand handled in `solve_demand.jl`.
=#

"""
    solve_consumption(
        archetype::ArchetypeBuilding,
        hvac_demand::Dict{Object,T} where {T<:SpineDataType};
        mod::Module = @__MODULE__,
    )

Solve the consumption of `abstract_process` included in the `archetype`,
based on the `hvac_demand`.

NOTE! The `mod` keyword changes from which Module data is accessed from,
`@__MODULE__` by default.

Currently produces only extremely simple estimates,
where each `abstract_process` is assumed to handle the demand on their output
nodes in full, regardless of what the other `abstract_processes` are doing.
Thus, the results cannot be used as is for systems with complimentary HVAC
processess, and the merit-order assumptions are left up to the user.

Heating and cooling are handled separately, though.
"""
function solve_consumption(
    archetype::ArchetypeBuilding,
    hvac_demand::Dict{Object,T} where {T<:SpineDataType};
    mod::Module = @__MODULE__,
)
    # Initialize a Dict for storing the results, and convert `hvac_demand`
    # into a parameter value for easier access.
    hvac_consumption = Dict()
    for (process, abstract_process) in archetype.abstract_processes
        # Figure out which node the process handles using `maximum_flows`
        output_nodes =
            getindex.(
                filter(
                    tup -> tup[1] == mod.direction(:to_node),
                    keys(abstract_process.maximum_flows),
                ),
                2,
            )

        # Calculate the consumption of the process.
        cons = reduce(
            +,
            hvac_demand[node] / abstract_process.coefficient_of_performance for
            node in output_nodes;
            init = 0.0,
        )

        # Remove negative consumption, as that indicates heating/cooling equipment
        # being used for cooling/heating inappropriately.
        hvac_consumption[process] = timedata_operation(v -> max(v, 0.0), cons)
    end
    return hvac_consumption
end
