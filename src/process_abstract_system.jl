#=
    process_abstract_system.jl

Contains functions for calculating the properties of abstract HVAC systems,
in preparation of conversion to large-scale energy system model input.
=#

"""
    process_abstract_system(process::BuildingProcessData; mod::Module = @__MODULE__)

Process a [`BuildingProcessData`](@ref) into an [`AbstractProcess`](@ref).

NOTE! The `mod` keyword changes from which Module data is accessed from,
`@__MODULE__` by default.

Essentially, combines the properties of a [`BuildingProcessData`](@ref)
into the bare minimum parameters required to describe a "process"
in our large-scale energy system modelling frameworks.
Performs the following steps:
1. COP mode indicated using COP sign, positive for heating and negative for cooling.
2. Scale COP to account for boundary between buildings and system link nodes. W -> MW conversion and accounting for the total number of processes.
3. Factor COP mode into the `maximum_flows_W` dictionary.
4. Return the components for [`AbstractProcess`](@ref).
"""
function process_abstract_system(process::BuildingProcessData; mod::Module = @__MODULE__)
    # COP sign used to indicate heating (positive) or cooling (negative).
    if process.coefficient_of_performance_mode == :heating
        coefficient_of_performance = process.coefficient_of_performance
    elseif process.coefficient_of_performance_mode == :cooling
        coefficient_of_performance = -process.coefficient_of_performance
    else
        @error """
        Unrecognized `coefficient_of_performance_mode` for `building_process` `$(process)`!
        Only `heating` and `cooling` supported!
        """
    end

    # Check if sys_link_node is in input/output nodes for the process and scale COP accordingly.
    if any( # Takes input from the system link node.
        in.(
            mod.building_process__direction__building_node(
                building_process = process.building_process,
                direction = mod.direction(:from_node),
            ),
            Ref(process.system_link_nodes),
        ),
    )
        coefficient_of_performance *= 1e6 / process.number_of_processes
    elseif any( # Produces output to the system link node.
        in.(
            mod.building_process__direction__building_node(
                building_process = process.building_process,
                direction = mod.direction(:to_node),
            ),
            Ref(process.system_link_nodes),
        ),
    )
        coefficient_of_performance *= process.number_of_processes / 1e6
    end

    # Maximum flows scaled using number of processes, sign used to indicate heating (positive) or cooling (negative).
    maximum_flows_W = Dict(
        (dir, node) =>
            ( # Cooling processes indicated using negative flows.
                dir == mod.direction(:to_node) &&
                process.coefficient_of_performance_mode == :cooling ? -1 : 1
            ) *
            ( # Scaling W -> MW
                !in(node, process.system_link_nodes) ? 1.0 :
                process.number_of_processes / 1e6
            ) *
            (
                process.maximum_power_base_W[(dir, node)] +
                process.maximum_power_gfa_scaled_W[(dir, node)]
            ) for (dir, node) in mod.building_process__direction__building_node(
            building_process = process.building_process,
        )
    )
    filter!(pair -> pair[2] != 0, maximum_flows_W)

    # Return the components for `AbstractProcess`.
    return process.number_of_processes, coefficient_of_performance, maximum_flows_W
end
