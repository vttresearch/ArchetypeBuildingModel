#=
    process_abstract_node.jl

Contains functions for processing the properties of
abstract lumped-capacitance thermal nodes, in preparation for conversion to
large-scale energy system model input.
=#

"""
    create_abstract_node_network(
        archetype::Object,
        scope::ScopeData,
        envelope::EnvelopeData,
        building_node_network::BuildingNodeNetwork,
        weather::WeatherData;
        mod::Module=@__MODULE__
    )

Process a `BuildingNodeNetwork` into an `AbstractNodeNetwork`.

The main purpose of this step is to account for weather-dependencies
contained in the `weather::WeatherData` input. However, the
`AbstractNodeNetwork` is also computationally useful,
e.g. for creating model-agnostic input for multiple
large-scale energy system models.
"""
function create_abstract_node_network(
    archetype::Object,
    scope::ScopeData,
    envelope::EnvelopeData,
    building_node_network::BuildingNodeNetwork,
    weather::WeatherData;
    mod::Module=@__MODULE__
)
    Dict(
        node => AbstractNode(
            archetype,
            scope,
            envelope,
            building_node_network,
            weather,
            node;
            mod=mod
        ) for
        node in keys(building_node_network)
    )::AbstractNodeNetwork
end


"""
    process_abstract_node(
        archetype::Object,
        scope::ScopeData,
        envelope::EnvelopeData,
        building_node_network::BuildingNodeNetwork,
        weather::WeatherData,
        node::Object;
        mod::Module=@__MODULE__
    )

Calculate the properties of an [`AbstractNode`](@ref) corresponding to the `node` in the `building_node_network`.

Combines all the individual parameters in [`BuildingNodeData`](@ref)s in [`BuildingNodeNetwork`](@ref)
into the bare minimum parameters required for modelling lumped-capacitance thermal nodes
in our energy system modelling frameworks.

Essentially, this function performs the following steps:
1. Sum all the thermal mass components together and convert into [Wh/K].
2. Sum all the self-discharge and ambient heat transfer components together.
3. Collect heat transfer coefficients between the interior air node and this one.
4. Update heat transfer coefficients based on user-defined coefficients.
5. Sum together the internal heat gains, solar gains, radiative sky heat losses, DHW demand, as well as the impact of ambient temperatures.
6. Return the components required for constructing an [`AbstractNode`](@ref).

**NOTE!** The ambient temperatures are accounted for via a combination of `self_discharge_coefficient_W_K`
and `external_load_W`, instead of  `heat_transfer_coefficients_W_K` on any ambient temperature nodes,
as illustrated by the equations below.
Typically, the heat balance equation in simplified lumped-capacitance thermal
models is cast as
```math
C_n \\frac{dT_n(t)}{dt} = H_{amb,n} \\left( T_{amb}(t) - T_{n}(t) \\right) + \\sum_{m \\in \\mathbb{N}} \\left[ H_{n,m} \\left( T_{m}(t) - T_{n}(t) \\right) \\right] + \\Sigma P_{n}(t) + \\Sigma \\Phi_{n}(t),
```
where `C_n` is the effective thermal mass of node `n`,
`T_{n}(t)` is the temperature of node `n` on time step `t`,
`T_{amb}(t)` is the ambient temperature,
`H_{amb,n}` is the conductance between ambient temperature and the node temperature,
`N` is the set of temperature nodes connected to node `n`,
`H_{n,m}` is the conductance between nodes `n` and `m`,
`∑P_{n}(t)` is the total impact of HVAC equipment,
and `∑Φ_{n}(t)` is the total effect of internal and solar heat gains.
However, large-scale energy system models rarely support ambient temperature `T_{amb}(t)`
as input data directly, requiring the above equation to be cast as
```math
C_n \\frac{dT_n(t)}{dt} = - H_{amb,n} T_{n}(t) + \\sum_{m \\in N} \\left[ H_{n,m} \\left( T_{m}(t) - T_{n}(t) \\right) \\right] + \\Sigma P_{n}(t) + \\left( H_{amb,n} T_{amb}(t) + \\Sigma \\Phi_{n}(t) \\right).
```
Now the `- H_{amb,n} T_{n}(t)` term can be interpreted as self-discharge losses,
while the `H_{amb,n} T_{amb}(t)` term can be bundled together with other external
influences on the node, both supported by typical large-scale energy system models.
Unfortunately, this has the side-effect of making the energy-system-model-level
input data quite unintuitive, but avoids the need to implement ambient-temperature-dependent
interactions in complicated energy system modelling frameworks.

**NOTE!** All heat transfer coefficients are forced to be symmetrical!
**NOTE!** Currently, radiative internal and solar gains are lost through the windows
in the building envelope.
"""
function process_abstract_node(
    archetype::Object,
    scope::ScopeData,
    envelope::EnvelopeData,
    building_node_network::BuildingNodeNetwork,
    weather::WeatherData,
    node::Object;
    mod::Module=@__MODULE__
)
    # Convenience access to the `BuildingNodeData` and the other nodes, making sure not to alter the original network.
    building_node_network = copy(building_node_network)
    node_data = pop!(building_node_network, node)

    # Total thermal mass of the node, in Wh/K for better scaling in energy system models.
    thermal_mass_Wh_K =
        (
            node_data.thermal_mass_base_J_K +
            node_data.thermal_mass_gfa_scaled_J_K +
            node_data.thermal_mass_interior_air_and_furniture_J_K +
            node_data.thermal_mass_structures_J_K
        ) / 3600

    # Total self-discharge coefficient from the node, accounting for ambient heat transfer.
    self_discharge_coefficient_W_K =
        node_data.self_discharge_base_W_K +
        node_data.self_discharge_gfa_scaled_W_K +
        node_data.heat_transfer_coefficient_structures_exterior_W_K +
        node_data.heat_transfer_coefficient_structures_ground_W_K +
        node_data.heat_transfer_coefficient_windows_W_K +
        node_data.heat_transfer_coefficient_ventilation_and_infiltration_W_K +
        node_data.heat_transfer_coefficient_thermal_bridges_W_K

    # Heat transfer coefficients from this node to connected nodes.
    heat_transfer_coefficients_W_K = mergewith(
        +,
        Dict( # First, heat transfer to interior air.
            n =>
                node_data.heat_transfer_coefficient_structures_interior_W_K *
                n_data.interior_air_and_furniture_weight for
            (n, n_data) in building_node_network
        ),
        Dict( # Force symmetry to the interior air heat transfer.
            n =>
                n_data.heat_transfer_coefficient_structures_interior_W_K *
                node_data.interior_air_and_furniture_weight for
            (n, n_data) in building_node_network
        ),
        node_data.heat_transfer_coefficients_base_W_K, # User-defined heat transfers
        node_data.heat_transfer_coefficients_gfa_scaled_W_K, # User-defined heat transfers
        Dict( # Force symmetry on user-defined heat transfer
            n => (
                get(n_data.heat_transfer_coefficients_base_W_K, node, 0.0) +
                get(n_data.heat_transfer_coefficients_gfa_scaled_W_K, node, 0.0)
            ) for (n, n_data) in building_node_network
        )
    )
    # And filter out zero heat transfer coefficients.
    filter!(pair -> pair[2] != 0, heat_transfer_coefficients_W_K)

    # Calculate the necessary solar heat gains
    solar_heat_gains_air_W = calculate_convective_solar_gains(
        archetype,
        scope,
        envelope,
        weather,
        node_data.interior_air_and_furniture_weight;
        mod=mod
    )
    solar_heat_gains_structures_W = calculate_radiative_solar_gains(
        archetype,
        node,
        scope,
        envelope,
        weather;
        mod=mod
    )
    solar_heat_gains_envelope_W = calculate_total_envelope_solar_gains(
        archetype,
        node,
        scope,
        envelope,
        weather;
        mod=mod
    )
    radiative_envelope_sky_losses_W = calculate_total_envelope_radiative_sky_losses(
        archetype,
        node,
        scope,
        envelope;
        mod=mod
    )

    # External load accounting for heat transfer with ambient conditions.
    external_load_W =
        (
            node_data.heat_transfer_coefficient_structures_exterior_W_K +
            node_data.heat_transfer_coefficient_windows_W_K +
            node_data.heat_transfer_coefficient_ventilation_and_infiltration_W_K +
            node_data.heat_transfer_coefficient_thermal_bridges_W_K
        ) * weather.ambient_temperature_K +
        node_data.heat_transfer_coefficient_structures_ground_W_K *
        weather.ground_temperature_K +
        node_data.internal_heat_gains_air_W +
        node_data.internal_heat_gains_structures_W +
        solar_heat_gains_air_W +
        solar_heat_gains_structures_W +
        solar_heat_gains_envelope_W -
        radiative_envelope_sky_losses_W -
        node_data.domestic_hot_water_demand_W

    # Return the properties of interest in the correct order for `AbstractNode`.
    return node_data.heating_set_point_K,
    node_data.cooling_set_point_K,
    node_data.maximum_temperature_deviation_K,
    node_data.minimum_temperature_deviation_K,
    thermal_mass_Wh_K,
    self_discharge_coefficient_W_K,
    heat_transfer_coefficients_W_K,
    external_load_W,
    node_data.interior_air_and_furniture_weight
end


"""
    calculate_window_solar_gains(
        archetype::Object,
        scope::ScopeData,
        envelope::EnvelopeData,
        weather::WeatherData;
        mod::Module = @__MODULE__,
    )

TODO: Revise documentation!

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
function calculate_window_solar_gains(
    archetype::Object,
    scope::ScopeData,
    envelope::EnvelopeData,
    weather::WeatherData;
    mod::Module=@__MODULE__
)
    mod.window_non_perpendicularity_correction_factor(building_archetype=archetype) *
    scope.total_normal_solar_energy_transmittance *
    envelope.window.surface_area_m2 *
    sum(
        mod.window_area_distribution_towards_cardinal_directions(
            building_archetype=archetype;
            cardinal_direction=dir
        ) * weather.total_effective_solar_irradiation_W_m2[sol_dir_map[dir]]
        for dir in old_solar_directions
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

TODO: Revise documentation!

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
\\Phi_\\text{sol,st} = R_\\text{e,st} U_\\text{ext,st} A_\\text{st} a_\\text{sol} \\left( I_\\text{diff} + F_\\text{shading} \\frac{\\sum_{d \\in D_\\text{st}} I_\\text{dir,d}}{\\sum_{d \\in D_\\text{st}} 1} \\right)
```
where `R_e,st` is the [exterior\\_resistance\\_m2K\\_W](@ref) of structure `st`,
`U_ext,st` is the [external\\_U\\_value\\_to\\_ambient\\_air\\_W\\_m2K](@ref) of structure `st`,
`A_st` is the surface area of the corresponding envelope structure (See [`EnvelopeData`](@ref)),
`a_sol` is the assumed [average\\_structural\\_solar\\_absorption\\_coefficient](@ref),
`I_diff` is the [diffuse\\_solar\\_irradiation\\_W\\_m2](@ref),
`F_shading` is the assumed [external\\_shading\\_coefficient](@ref),
`d` represents either horizontal or cardinal directions,
and `I_dir,d` is the [direct\\_solar\\_irradiation\\_W\\_m2](@ref).

NOTE! The walls are assumed to be distributed equally towards all the cardinal directions.
"""
function calculate_envelope_solar_gains(
    archetype::Object,
    scope::ScopeData,
    envelope::EnvelopeData,
    weather::WeatherData;
    mod::Module=@__MODULE__
)
    # Define applicable structure types and directions.
    st_dirs = Dict(
        mod.structure_type(:exterior_wall) => :vertical,
        mod.structure_type(:light_exterior_wall) => :vertical,
        mod.structure_type(:roof) => :horizontal,
    )

    # Calculate the net solar gains and sky losses for each structure type.
    return Dict(
        st => (
            mod.exterior_resistance_m2K_W(structure_type=st) *
            scope.structure_data[st].external_U_value_to_ambient_air_W_m2K *
            getfield(envelope, st.name).surface_area_m2 *
            mod.average_structural_solar_absorption_coefficient(
                building_archetype=archetype,
            ) *
            weather.total_effective_solar_irradiation_W_m2[dir]
        ) for (st, dir) in st_dirs
    )
end


"""
    calculate_envelope_radiative_sky_losses(
        archetype::Object,
        scope::ScopeData,
        envelope::EnvelopeData;
        mod::Module = @__MODULE__,
    )

TODO: REVISE DOCUMENTATION!

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
    mod::Module=@__MODULE__
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
            mod.exterior_resistance_m2K_W(structure_type=st) *
            scope.structure_data[st].external_U_value_to_ambient_air_W_m2K *
            getfield(envelope, st.name).surface_area_m2 *
            F_sky *
            mod.external_radiative_surface_heat_transfer_coefficient_W_m2K(
                building_archetype=archetype,
            ) *
            mod.average_apparent_sky_temperature_difference_K(
                building_archetype=archetype,
            )
        ) for (st, F_sky) in st_F_sky
    )
end



"""
    calculate_convective_solar_gains(
        archetype::Object
        loads::LoadsData,
        interior_weight::Real;
        mod::Module = @__MODULE__,
    )

TODO: Revise documentation!

Calculate the convective solar heat gains through windows on the `node` in [W].

NOTE! The `mod` keyword changes from which Module data is accessed from,
`@__MODULE__` by default.

Essentially, takes the given solar heat gain profile in `loads` and
multiplies it with the share of interior air on this `node` as well as the
assumed convective fraction of solar heat gains.
```math
\\Phi_\\text{sol,conv,n} = w_\\text{int,n} f_\\text{sol,conv} \\Phi_\\text{sol}
```
where `w_int,n` is the [interior\\_air\\_and\\_furniture\\_weight](@ref) of this node,
`f_sol,conv` is the assumed [solar\\_heat\\_gain\\_convective\\_fraction](@ref),
and `Φ_sol` are the total solar heat gains into the building.
See [`calculate_total_solar_gains`](@ref) for how the total solar heat gains
are calculated.
"""
function calculate_convective_solar_gains(
    archetype::Object,
    scope::ScopeData,
    envelope::EnvelopeData,
    weather::WeatherData,
    interior_weight::Real;
    mod::Module=@__MODULE__
)
    calculate_window_solar_gains(
        archetype,
        scope,
        envelope,
        weather;
        mod=mod
    ) *
    interior_weight *
    mod.solar_heat_gain_convective_fraction(building_archetype=archetype)
end


"""
    calculate_radiative_solar_gains(
        archetype::Object,
        node::Object,
        envelope::EnvelopeData,
        loads::LoadsData,
        total_structure_area_m2::Real;
        mod::Module = @__MODULE__,
    )

TODO: Revise documentation

Calculate the radiative solar heat gains through windows on the `node` in [W].

NOTE! The `mod` keyword changes from which Module data is accessed from,
`@__MODULE__` by default.

Essentially, takes the given solar heat gain profile in `loads`
and multiplies it with the assumed radiative fraction of solar heat gains.
The radiative heat gains are assumed to be distributed across the structures 
simply based on their relative surface areas.
Note that currently, radiative solar heat gains are partially lost through windows!
```math
\\Phi_\\text{sol,rad,n} = (1 - f_\\text{sol,conv}) \\frac{\\sum_{\\text{st} \\in n} w_\\text{n,st} A_\\text{st}}{\\sum_{\\text{st}} A_\\text{st}} \\Phi_\\text{sol}
```
where `f_sol,conv` is the assumed [solar\\_heat\\_gain\\_convective\\_fraction](@ref),
`st` is the [structure\\_type](@ref) and `n` is this [building\\_node](@ref),
`w_n,st` is the [structure\\_type\\_weight](@ref) of the structure `st` on this node,
`A_st` is the surface area of structure `st`,
and `Φ_sol` are the total solar heat gains into the building.
See [`calculate_total_solar_gains`](@ref) for how the total solar heat gains
are calculated.
"""
function calculate_radiative_solar_gains(
    archetype::Object,
    node::Object,
    scope::ScopeData,
    envelope::EnvelopeData,
    weather::WeatherData;
    mod::Module=@__MODULE__
)
    calculate_window_solar_gains(
        archetype,
        scope,
        envelope,
        weather;
        mod=mod
    ) *
    (1 - mod.solar_heat_gain_convective_fraction(building_archetype=archetype)) *
    sum(
        mod.structure_type_weight(building_node=node, structure_type=st) *
        getfield(envelope, st.name).surface_area_m2 / envelope.total_structure_area_m2 for
        st in mod.building_node__structure_type(building_node=node);
        init=0
    )
end


"""
    calculate_total_envelope_solar_gains(
        node::Object,
        loads::LoadsData;
        mod::Module = @__MODULE__,
    )

TODO: Revise documentation!

Calculate the total solar heat gain [W] through the opaque envelope on this node.

NOTE! The `mod` keyword changes from which Module data is accessed from,
`@__MODULE__` by default.

```math
\\Phi_\\text{env,n} = \\sum_{st \\in n} w_\\text{n,st} \\Phi_\\text{sol,st},
```
`w_n,st` is the [structure\\_type\\_weight](@ref) of the structure `st` on this node `n`,
and `Φ_sol,st` are the heat gains through envelope structures using [`calculate_envelope_solar_gains`](@ref).
"""
function calculate_total_envelope_solar_gains(
    archetype::Object,
    node::Object,
    scope::ScopeData,
    envelope::EnvelopeData,
    weather::WeatherData;
    mod::Module=@__MODULE__
)
    envelope_solar_gains_W = calculate_envelope_solar_gains(
        archetype, scope, envelope, weather; mod=mod
    )
    sum(
        mod.structure_type_weight(building_node=node, structure_type=st) *
        get(envelope_solar_gains_W, st, 0.0) for
        st in mod.building_node__structure_type(building_node=node);
        init=0.0
    )
end


"""
    calculate_total_envelope_radiative_sky_losses(
        node::Object,
        loads::LoadsData;
        mod::Module = @__MODULE__,
    )

Calculate the total radiative envelope sky heat losses [W] for this node.

TODO: Revise documentation!

NOTE! The `mod` keyword changes from which Module data is accessed from,
`@__MODULE__` by default.

```math
\\Phi_\\text{sky,n} = \\sum_{st \\in n} w_\\text{n,st} \\Phi_\\text{sky,st},
```
`w_n,st` is the [structure\\_type\\_weight](@ref) of the structure `st` on this node `n`,
and `Φ_sky,st` are the envelope radiative sky heat losses using [`calculate_envelope_radiative_sky_losses`](@ref).
"""
function calculate_total_envelope_radiative_sky_losses(
    archetype::Object,
    node::Object,
    scope::ScopeData,
    envelope::EnvelopeData;
    mod::Module=@__MODULE__
)
    envelope_radiative_sky_losses_W = calculate_envelope_radiative_sky_losses(
        archetype, scope, envelope; mod=mod
    )
    sum(
        mod.structure_type_weight(building_node=node, structure_type=st) *
        get(envelope_radiative_sky_losses_W, st, 0.0) for
        st in mod.building_node__structure_type(building_node=node);
        init=0.0
    )
end