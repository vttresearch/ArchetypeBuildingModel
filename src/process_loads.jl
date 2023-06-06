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
4. Calculates total solar gains through windows using [`calculate_total_solar_gains`](@ref).
5. Calculates solar gains through envelope structures using [`calculate_envelope_solar_gains`](@ref).
6. Calculates envelope radiative sky losses using [`calculate_envelope_radiative_sky_losses`](@ref).
7. Returns the calculated DHW demand, internal gains, solar gains for windows and the envelope, and envelope radiative sky losses.
"""
function process_building_loads(
    archetype::Object,
    scope::ScopeData,
    envelope::EnvelopeData,
    weather::WeatherData;
    mod::Module = @__MODULE__
)
    # Find the `building_loads` connected to the `archetype`
    loads = first(mod.building_archetype__building_loads(building_archetype = archetype))

    # Calculate loads
    dhw_demand = calculate_total_dhw_demand(loads, scope; mod = mod)
    internal_gains = calculate_total_internal_heat_loads(loads, scope; mod = mod)
    solar_gains =
        calculate_total_solar_gains(archetype, scope, envelope, weather; mod = mod)
    envelope_solar_gains =
        calculate_envelope_solar_gains(archetype, scope, envelope, weather; mod = mod)
    envelope_sky_losses =
        calculate_envelope_radiative_sky_losses(archetype, scope, envelope; mod = mod)

    return dhw_demand,
    internal_gains,
    solar_gains,
    envelope_solar_gains,
    envelope_sky_losses
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
    mod::Module = @__MODULE__
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
    mod::Module = @__MODULE__
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
Solar heat gains through the envelope are handled via [`calculate_envelope_solar_gains`](@ref).

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
    mod::Module = @__MODULE__
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
            ) * weather.direct_solar_irradiation_W_m2[dir] for dir in solar_directions
        )
    )
end


"""
    calculate_envelope_solar_gains(
        archetype::Object,
        scope::ScopeData,
        envelope::EnvelopeData,
        weather::WeatherData;
        mod::Module = @__MODULE__,
    )

Calculate the solar heat gains [W] per envelope [structure\\_type](@ref).

NOTE! The `mod` keyword changes from which Module data is accessed from,
`@__MODULE__` by default.

Loosely based on EN ISO 52016-1:2017 Section 6.5.6.3.5,
approximately accounts for the incident solar irradiation, while
approximate radiative sky heat losses from the envelope are handled via
[`calculate_envelope_radiative_sky_losses`](@ref).
Solar heat gains through windows are handled via [`calculate_total_solar_gains`](@ref).

Essentially, the exterior of the surface is assumed to have negligible thermal mass,
and the impact of the incident irradiation is applied to the structural node directly.
The solar gains `Φ_sol,st` need to be calculated separately for
each [structure\\_type](@ref), as they are dependent on the exterior surface resistance:
```math
\\Phi_\\text{sol,st} = R_\\text{e,st} U_\\text{ext,st} A_\\text{st} a_\\text{sol} \\left( I_\\text{diff} + F_\\text{shading} \\sum_{d \\in D_\\text{st}} I_\\text{dir,d} \\right)
```
where `R_e,st` is the [exterior\\_resistance\\_m2K\\_W](@ref) of structure `st`,
`U_ext,st` is the [external\\_U\\_value\\_to\\_ambient\\_air\\_W\\_m2K](@ref) of structure `st`,
`A_st` is the surface area of the corresponding envelope structure (See [`EnvelopeData`](@ref)),
`a_sol` is the assumed [average\\_structural\\_solar\\_absorption\\_coefficient](@ref),
`I_diff` is the [diffuse\\_solar\\_irradiation\\_W\\_m2](@ref),
`F_shading` is the assumed [external\\_shading\\_coefficient](@ref),
`d` represents either horizontal or cardinal directions,
and `I_dir,d` is the [direct\\_solar\\_irradiation\\_W\\_m2](@ref).
"""
function calculate_envelope_solar_gains(
    archetype::Object,
    scope::ScopeData,
    envelope::EnvelopeData,
    weather::WeatherData;
    mod::Module = @__MODULE__
)
    # Define applicable structure types and directions.
    st_dirs = Dict(
        mod.structure_type(:exterior_wall) => [:north, :east, :south, :west],
        mod.structure_type(:light_exterior_wall) => [:north, :east, :south, :west],
        mod.structure_type(:roof) => [:horizontal],
    )

    # Calculate the net solar gains and sky losses for each structure type.
    return Dict(
        st => (
            mod.exterior_resistance_m2K_W(structure_type = st) *
            scope.structure_data[st].external_U_value_to_ambient_air_W_m2K *
            getfield(envelope, st.name).surface_area_m2 *
            mod.average_structural_solar_absorption_coefficient(
                building_archetype = archetype,
            ) *
            (
                weather.diffuse_solar_irradiation_W_m2 +
                mod.external_shading_coefficient(building_archetype = archetype) *
                sum(weather.direct_solar_irradiation_W_m2[dir] for dir in dirs)
            )
        ) for (st, dirs) in st_dirs
    )
end


"""
    calculate_envelope_radiative_sky_losses(
        archetype::Object,
        scope::ScopeData,
        envelope::EnvelopeData;
        mod::Module = @__MODULE__,
    )

Calculate the envelope radiative sky losses [W] per [structure\\_type](@ref).

NOTE! The `mod` keyword changes from which Module data is accessed from,
`@__MODULE__` by default.

Loosely based on EN ISO 52016-1:2017 Section 6.5.13.3,
approximately accounts for the radiative sky heat losses of the building envelope.
Solar heat gains through the building envelope are handled via
[`calculate_envelope_solar_gains`](@ref).

Essentially, the exterior of the surface is assumed to have negligible thermal mass,
and the impact of radiative sky losses are applied to the structural node directly.
The radiative sky heat losses `Φ_sky,st` are calculated for each structure as
```math
\\Phi_\\text{sky,st} = R_\\text{e,st} U_\\text{ext,st} A_\\text{st} F_\\text{sky,st} h_\\text{re} \\Delta T_\\text{sky}
```
where `R_e,st` is the [exterior\\_resistance\\_m2K\\_W](@ref) of structure `st`,
`U_ext,st` is the [external\\_U\\_value\\_to\\_ambient\\_air\\_W\\_m2K](@ref) of structure `st`,
`A_st` is the surface area of the corresponding envelope structure (See [`EnvelopeData`](@ref)),
`F_sky,st` is the assumed sky view factor *(hardcoded for now)* for structure `st`,
`h_re` is the assumed [external\\_radiative\\_surface\\_heat\\_transfer\\_coefficient\\_W\\_m2K](@ref),
and `ΔT_sky` is the assumed [average\\_apparent\\_sky\\_temperature\\_difference\\_K](@ref).
"""
function calculate_envelope_radiative_sky_losses(
    archetype::Object,
    scope::ScopeData,
    envelope::EnvelopeData;
    mod::Module = @__MODULE__
)
    # Define envelope structure types and their view factors to the sky, based on EN ISO 52016-1:2017 Table B.18.
    st_F_sky = Dict(
        mod.structure_type(:exterior_wall) => 0.5,
        mod.structure_type(:light_exterior_wall) => 0.5,
        mod.structure_type(:roof) => 1.0,
    )

    # Calculate the radiative sky heat losses for each structure
    return Dict(
        st => (
            mod.exterior_resistance_m2K_W(structure_type = st) *
            scope.structure_data[st].external_U_value_to_ambient_air_W_m2K *
            getfield(envelope, st.name).surface_area_m2 *
            F_sky *
            mod.external_radiative_surface_heat_transfer_coefficient_W_m2K(
                building_archetype = archetype,
            ) *
            mod.average_apparent_sky_temperature_difference_K(
                building_archetype = archetype,
            )
        ) for (st, F_sky) in st_F_sky
    )
end
