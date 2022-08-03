#=
    process_loads.jl

Contains functions for processing the domestic hot water demand, internal heat gains,
and solar gains for the archetype buildings.
=#


"""
    process_building_loads(
        archetype::Object,
        scope::ScopeData,
        envelope::EnvelopeData,
        weather::WeatherData;
        mod::Module = @__MODULE__,
    )

Calculate the domestic hot water demand, internal and solar heat gains for the archetype building.

NOTE! The `mod` keyword changes from which Module data is accessed from,
`@__MODULE__` by default.

Essentially, performs the following steps:
1. Finds the `building_loads` object corresponding to the `building_archetype`.
2. Calculates total domestic hot water (DHW) demand using [`calculate_total_dhw_demand`](@ref).
3. Calculates total internal heat loads using [`calculate_total_internal_heat_loads`](@ref).
4. Calculates total solar gains using [`calculate_total_solar_gains`](@ref).
5. Returns the calculated DHW demand, internal gains, and solar gains.
"""
function process_building_loads(
    archetype::Object,
    scope::ScopeData,
    envelope::EnvelopeData,
    weather::WeatherData;
    mod::Module = @__MODULE__,
)
    # Find the `building_loads` connected to the `archetype`
    loads = first(mod.building_archetype__building_loads(building_archetype = archetype))

    # Calculate loads
    dhw_demand = calculate_total_dhw_demand(loads, scope; mod = mod)
    internal_gains = calculate_total_internal_heat_loads(loads, scope; mod = mod)
    solar_gains =
        calculate_total_solar_gains(archetype, scope, envelope, weather; mod = mod)

    return dhw_demand, internal_gains, solar_gains
end


"""
    calculate_total_dhw_demand(loads::Object, scope::ScopeData; mod::Module = @__MODULE__)

Calculate the total domestic hot water (DHW) demand.

NOTE! The `mod` keyword changes from which Module data is accessed from,
`@__MODULE__` by default.

Essentially, simply adds the fixed base DHW demand and
the gross-floor-area-scaling DHW demand together
for the given `loads` and `scope`.
```math
\\Phi_\\text{DHW} = \\Phi_\\text{DHW,base} + \\phi_\\text{DHW,gfa} A_\\text{gfa}
```
where `Φ_DHW,base` is the assumed [domestic\\_hot\\_water\\_demand\\_base\\_W](@ref),
`Φ_DHW,gfa` is the assumed [domestic\\_hot\\_water\\_demand\\_gfa\\_scaling\\_W\\_m2](@ref),
and `A_gfa` is the gross-floor area of the building.
"""
function calculate_total_dhw_demand(
    loads::Object,
    scope::ScopeData;
    mod::Module = @__MODULE__,
)
    mod.domestic_hot_water_demand_base_W(building_loads = loads) +
    scope.average_gross_floor_area_m2_per_building *
    mod.domestic_hot_water_demand_gfa_scaling_W_m2(building_loads = loads)
end


"""
    calculate_total_internal_heat_loads(
        loads::Object,
        scope::ScopeData;
        mod::Module = @__MODULE__
    )

Calculate the total internal heat gains.

NOTE! The `mod` keyword changes from which Module data is accessed from,
`@__MODULE__` by default.

Essentially, simply adds the fixed base internal gains and
the gross-floor-area-scaling DHW demand together
for the given `loads` and `scope`.
```math
\\Phi_\\text{int} = \\Phi_\\text{int,base} + \\phi_\\text{int,gfa} A_\\text{gfa}
```
where `Φ_int,base` is the assumed [internal\\_heat\\_loads\\_base\\_W](@ref),
`Φ_int,gfa` is the assumed [internal\\_heat\\_loads\\_gfa\\_scaling\\_W\\_m2](@ref),
and `A_gfa` is the gross-floor area of the building.
"""
function calculate_total_internal_heat_loads(
    loads::Object,
    scope::ScopeData;
    mod::Module = @__MODULE__,
)
    mod.internal_heat_loads_base_W(building_loads = loads) +
    scope.average_gross_floor_area_m2_per_building *
    mod.internal_heat_loads_gfa_scaling_W_m2(building_loads = loads)
end


"""
    calculate_total_solar_gains(
        archetype::Object,
        scope::ScopeData,
        envelope::EnvelopeData,
        weather::WeatherData;
        mod::Module = @__MODULE__,
    )

Calculate the total solar heat gains through the windows.

NOTE! The `mod` keyword changes from which Module data is accessed from,
`@__MODULE__` by default.

Loosely based on EN ISO 52016-1:2017 6.5.13.2,
accounting for the solar heat gains through the windows of the building.
Solar heat gains through the opaque parts of the building envelope are neglegted.
The varying angle of incidence of the irradiation on the windows is accounted for
using a very simple average non-perpendicularity factor.
The frame-area fraction of the windows is accounted for in the solar energy transmittance,
and window area distribution is handled using the shares towards cardinal directions.
```math
\\Phi_\\text{sol} = f_\\text{np} g_\\text{gl} A_\\text{w} \\left( I_\\text{diff} + F_\\text{shading} \\sum_{d \\in \\text{N,E,S,W}} w_d I_\\text{dir,d} \\right)
```
where `f_np` is the assumed [window\\_non\\_perpendicularity\\_correction\\_factor](@ref),
`g_gl` is the [total\\_normal\\_solar\\_energy\\_transmittance](@ref) of the glazing,
`A_w` is the area of the windows,
`I_diff` is the [diffuse\\_solar\\_irradiation\\_W\\_m2](@ref),
`F_shading` is the assumed [external\\_shading\\_coefficient](@ref),
`d` represents the cardinal directions,
`w_d` is the [window\\_area\\_distribution\\_towards\\_cardinal\\_directions](@ref),
and `I_dir,d` is the [direct\\_solar\\_irradiation\\_W\\_m2](@ref).
"""
function calculate_total_solar_gains(
    archetype::Object,
    scope::ScopeData,
    envelope::EnvelopeData,
    weather::WeatherData;
    mod::Module = @__MODULE__,
)
    mod.window_non_perpendicularity_correction_factor(building_archetype = archetype) *
    scope.total_normal_solar_energy_transmittance *
    envelope.window.surface_area_m2 *
    (
        weather.diffuse_solar_irradiation_W_m2 +
        mod.external_shading_coefficient(building_archetype = archetype) * sum(
            mod.window_area_distribution_towards_cardinal_directions(
                building_archetype = archetype;
                cardinal_direction = dir,
            ) * weather.direct_solar_irradiation_W_m2[dir] for
            dir in [:north, :east, :south, :west]
        )
    )
end
