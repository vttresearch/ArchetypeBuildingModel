#=
    process_system.jl

Contains functions for calculating the properties of HVAC systems.
=#


"""
    process_building_system(
        archetype::Object,
        process::Object,
        scope::ScopeData,
        weather::WeatherData;
        mod::Module = @__MODULE__,
    )

Calculates the properties required for creating a [`BuildingProcessData`](@ref).

NOTE! The `mod` keyword changes from which Module data is accessed from,
`@__MODULE__` by default.

Essentially performs the following steps:
1. Fetch the `number_of_processes`, `system_link_nodes`, and `COP_mode` from the definitions.
2. Calculate the coefficient of performance using the [`calculate_cop`](@ref) function.
3. Fetch the defined maximum power flows.
4. Return the components for [`BuildingProcessData`](@ref).
"""
function process_building_system(
    archetype::Object,
    process::Object,
    scope::ScopeData,
    weather::WeatherData;
    mod::Module = @__MODULE__,
)
    # Record the number of processes and system link nodes for input/output scaling of AbstractProcess
    number_of_processes = scope.number_of_buildings
    system_link_nodes =
        mod.building_archetype__system_link_node(building_archetype = archetype)

    # Fetch/calculate the different elements of the coefficient of performance.
    COP_mode = mod.coefficient_of_performance_mode(building_process = process)
    COP = calculate_cop(weather, COP_mode, process; mod = mod)

    # Fetch and calculate the process flow maximum power terms.
    maximum_power_basis_W = Dict(
        (dir, node) => mod.maximum_power_base_W(
            building_process = process,
            direction = dir,
            building_node = node,
        ) for (dir, node) in
        mod.building_process__direction__building_node(building_process = process)
    )
    maximum_power_gfa_scaled_W = Dict(
        (dir, node) =>
            scope.average_gross_floor_area_m2_per_building *
            mod.maximum_power_gfa_scaling_W_m2(
                building_process = process,
                direction = dir,
                building_node = node,
            ) for (dir, node) in
        mod.building_process__direction__building_node(building_process = process)
    )

    # Return the components of `BuildingProcessData`.
    return system_link_nodes,
    COP,
    COP_mode,
    maximum_power_basis_W,
    maximum_power_gfa_scaled_W,
    number_of_processes
end


"""
    calculate_cop(
        weather::WeatherData,
        COP_mode::Symbol,
        process::Object;
        mod::Module = @__MODULE__,
    )

Calculate the potentially time-varying coefficient of performance.

NOTE! The `mod` keyword changes from which Module data is accessed from,
`@__MODULE__` by default.

Uses the given `COP_mode`, `sink_temperature_K`, and `source_temperature_K` to
calculate a potentially weather-dependent coefficient of performance.
Based on the exergetic approach in EN 15316-4-2:2017 Annex D,
where a known reference COP is extrapolated to unknown temperature ranges
using the Carnot COP.
```math
\\text{COP}_{T_\\text{in}, T_\\text{out}} = \\text{COP}_\\text{base} \\frac{T_\\text{out}}{T_\\text{out} - T_\\text{in}}, \\\\
\\text{where } \\text{COP}_\\text{base} = \\frac{\\text{COP}_\\text{ref}}{\\text{COP}_\\text{carnot,ref}}
```

If the source and sink temperatures are undefined, returns
the defined base COP as the COP is assumed to be weather independent.
See the `archetype_definitions.json` default values for `sink_temperature_K`
and `source_temperature_K` for the different options how the COP can be made
dependent on the `weather`.
"""
function calculate_cop(
    weather::WeatherData,
    COP_mode::Symbol,
    process::Object;
    mod::Module = @__MODULE__,
)
    # Check if source and sink temperatures are defined, return 1.0 if not.
    if (
        mod.coefficient_of_performance_sink_temperature_K(building_process = process) isa
        Map{Symbol,SpineInterface.NothingParameterValue} &&
        mod.coefficient_of_performance_source_temperature_K(building_process = process) isa
        Map{Symbol,SpineInterface.NothingParameterValue}
    )
        COP_multiplier = 1.0
    else
        # Initialize a Dict for sink and source temperatures.
        temperature_K = Dict{Symbol,SpineDataType}()

        # Convert into `MapParameterValue` for convenience.
        for (key, temp) in zip(
            [:sink, :source],
            [
                mod.coefficient_of_performance_sink_temperature_K,
                mod.coefficient_of_performance_source_temperature_K,
            ],
        )

            # Find the proper control temperature
            if temp(building_process = process, property = :temperature_K) == :ambient
                control_temperature_K = weather.ambient_temperature_K
            elseif temp(building_process = process, property = :temperature_K) == :ground
                control_temperature_K = weather.ground_temperature_K
            else
                control_temperature_K =
                    temp(building_process = process, property = :temperature_K)
            end

            # Determine the heating curve if possible
            temps = [
                temp(
                    building_process = process,
                    property = :heating_curve_control_temperature_min_K,
                ),
                temp(
                    building_process = process,
                    property = :heating_curve_control_temperature_max_K,
                ),
                temp(
                    building_process = process,
                    property = :heating_curve_output_temperature_min_K,
                ),
                temp(
                    building_process = process,
                    property = :heating_curve_output_temperature_max_K,
                ),
            ]
            if !all(isnothing.(temps))
                lin = LinearInterpolation(
                    [temps[1], temps[2]],
                    [temps[4], temps[3]];
                    extrapolation_bc = Flat(),
                )
                heating_curve = x -> lin(x) # Generic function for linear interpolation.
            else
                heating_curve = x -> x # Generic function to return whatever it receives.
            end

            # Calculate the temperature using the control temperature and the heating curve.
            if control_temperature_K isa Union{TimePattern,TimeSeries}
                temperature_K[key] =
                    timedata_operation(heating_curve, control_temperature_K)
            else
                temperature_K[key] = heating_curve(control_temperature_K)
            end
        end

        # Carnot COP calculated with assumed minimum temperature delta.
        denom = temperature_K[:sink] - temperature_K[:source]
        if denom isa Union{TimePattern,TimeSeries}
            denom = timedata_operation(
                max,
                denom,
                mod.coefficient_of_performance_minimum_temperature_delta(
                    building_process = process,
                ),
            )
        else
            denom = max(
                denom,
                mod.coefficient_of_performance_minimum_temperature_delta(
                    building_process = process,
                ),
            )
        end

        # Calculate the COP multiplier depending on `COP_mode`
        if COP_mode == :heating
            COP_multiplier = temperature_K[:sink] / denom
        elseif COP_mode == :cooling
            COP_multiplier = temperature_K[:source] / denom
        else
            @error """
            Unrecognized `coefficient_of_performance_mode`!
            Only `:heating` and `:cooling` are recognized modes!
            """
            COP_multiplier = nothing
        end
    end

    # Return the coefficient of performance.
    return COP_multiplier * mod.coefficient_of_performance_base(building_process = process)
end


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
    return coefficient_of_performance, maximum_flows_W
end
