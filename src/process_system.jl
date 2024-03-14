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
1. Fetch the `system_link_nodes`, and `COP_mode` from the definitions.
2. Calculate the coefficient of performance using the [`calculate_cop`](@ref) function.
3. Fetch the defined maximum power flows.
4. Return the components for [`BuildingProcessData`](@ref).
"""
function process_building_system(
    archetype::Object,
    process::Object,
    scope::ScopeData,
    weather::WeatherData;
    mod::Module=@__MODULE__
)
    # Fetch system link nodes
    system_link_nodes =
        mod.building_archetype__system_link_node(building_archetype=archetype)

    # Fetch/calculate the different elements of the coefficient of performance.
    COP_mode = mod.coefficient_of_performance_mode(building_process=process)
    COP = calculate_cop(weather, COP_mode, process; mod=mod)

    # Fetch and calculate the process flow maximum power terms.
    maximum_power_basis_W = Dict(
        (dir, node) => mod.maximum_power_base_W(
            building_process=process,
            direction=dir,
            building_node=node,
        ) for (dir, node) in
        mod.building_process__direction__building_node(building_process=process)
    )
    maximum_power_gfa_scaled_W = Dict(
        (dir, node) =>
            scope.average_gross_floor_area_m2_per_building *
            mod.maximum_power_gfa_scaling_W_m2(
                building_process=process,
                direction=dir,
                building_node=node,
            ) for (dir, node) in
        mod.building_process__direction__building_node(building_process=process)
    )
    # Calculate the total maximum flows for convenience.
    maximum_flows_W = mergewith(
        +,
        maximum_power_basis_W,
        maximum_power_gfa_scaled_W,
    )
    filter!(pair -> pair[2] != 0, maximum_flows_W)

    # Return the components of `BuildingProcessData`.
    return system_link_nodes,
    COP,
    COP_mode,
    maximum_power_basis_W,
    maximum_power_gfa_scaled_W,
    maximum_flows_W
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
    mod::Module=@__MODULE__
)
    # Check if source and sink temperatures are defined, return 1.0 if not.
    if (
        mod.coefficient_of_performance_sink_temperature_K(building_process=process) isa
        Map{Symbol,Nothing} &&
        mod.coefficient_of_performance_source_temperature_K(building_process=process) isa
        Map{Symbol,Nothing}
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
            if temp(building_process=process, property=:temperature_K) == :ambient
                control_temperature_K = weather.ambient_temperature_K
            elseif temp(building_process=process, property=:temperature_K) == :ground
                control_temperature_K = weather.ground_temperature_K
            else
                control_temperature_K =
                    temp(building_process=process, property=:temperature_K)
            end

            # Determine the heating curve if possible
            temps = [
                temp(
                    building_process=process,
                    property=:heating_curve_control_temperature_min_K,
                ),
                temp(
                    building_process=process,
                    property=:heating_curve_control_temperature_max_K,
                ),
                temp(
                    building_process=process,
                    property=:heating_curve_output_temperature_min_K,
                ),
                temp(
                    building_process=process,
                    property=:heating_curve_output_temperature_max_K,
                ),
            ]
            if !all(isnothing.(temps))
                lin = LinearInterpolation(
                    [temps[1], temps[2]],
                    [temps[4], temps[3]];
                    extrapolation_bc=Flat()
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
                    building_process=process,
                ),
            )
        else
            denom = max(
                denom,
                mod.coefficient_of_performance_minimum_temperature_delta(
                    building_process=process,
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
    return COP_multiplier * mod.coefficient_of_performance_base(building_process=process)
end
