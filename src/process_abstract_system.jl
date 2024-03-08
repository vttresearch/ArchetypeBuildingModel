#=
    process_abstract_system.jl

Contains functions for calculating the properties of abstract HVAC systems,
in preparation of conversion to large-scale energy system model input.
=#

"""
    process_abstract_system(process::BuildingProcessData; mod::Module = @__MODULE__)

Process a [`BuildingProcessData`](@ref) into an [`AbstractProcess`](@ref).

TODO: Revise documentation!

NOTE! The `mod` keyword changes from which Module data is accessed from,
`@__MODULE__` by default.

Essentially, combines the properties of a [`BuildingProcessData`](@ref)
into the bare minimum parameters required to describe a "process"
in our large-scale energy system modelling frameworks.
Performs the following steps:
1. COP mode indicated using COP sign, positive for heating and negative for cooling.
2. Scale COP to account for boundary between buildings and system link nodes. W -> MW conversion and accounting for the total number of processes.
3. Factor COP mode into the `maximum_flows` dictionary.
4. Return the components for [`AbstractProcess`](@ref).
"""
function process_abstract_system(process::BuildingProcessData; mod::Module=@__MODULE__)
    # Maximum flows scaled using number of processes, sign used to indicate heating (positive) or cooling (negative).
    maximum_flows = Dict(
        (dir, node) => (
            process.maximum_power_base_W[(dir, node)] +
            process.maximum_power_gfa_scaled_W[(dir, node)]
        ) for (dir, node) in mod.building_process__direction__building_node(
            building_process=process.building_process,
        )
    )
    filter!(pair -> pair[2] != 0, maximum_flows)

    # Return the components for `AbstractProcess`.
    return process.coefficient_of_performance_mode,
    process.coefficient_of_performance,
    maximum_flows
end
