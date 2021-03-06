#=
    process_node.jl

Contains functions for processing the properties of lumped-capacitance thermal nodes.
=#

"""
    create_building_node_network(
        archetype::Object,
        fabrics::Object,
        systems::Object,
        scope::ScopeData,
        envelope::EnvelopeData,
        loads::LoadsData;
        mod::Module = @__MODULE__,
    )

Map all `archetype` `building_nodes` to their [`BuildingNodeData`](@ref)s.

NOTE! The `mod` keyword changes from which Module data is accessed from,
`@__MODULE__` by default.

Essentially, loops over the `building_fabrics__building_node` and
`building_systems__building_node` relationships for the desired `archetype`,
and collects all the `building_node`s and [`BuildingNodeData`](@ref)s
into a [`BuildingNodeNetwork`](@ref) dictionary.
"""
function create_building_node_network(
    archetype::Object,
    fabrics::Object,
    systems::Object,
    scope::ScopeData,
    envelope::EnvelopeData,
    loads::LoadsData;
    mod::Module = @__MODULE__,
)
    merge(
        Dict{Object,BuildingNodeData}(
            node => BuildingNodeData(archetype, node, scope, envelope, loads; mod = mod) for
            node in mod.building_fabrics__building_node(building_fabrics = fabrics)
        ),
        Dict{Object,BuildingNodeData}(
            node => BuildingNodeData(archetype, node, scope, envelope, loads; mod = mod) for
            node in mod.building_systems__building_node(building_systems = systems)
        ),
    )::BuildingNodeNetwork
end


"""
    process_building_node(
        node::Object,
        scope::ScopeData,
        envelope::EnvelopeData,
        loads::LoadsData;
        mod::Module = @__MODULE__,
    )

Calculates the properties of a `node` using the provided `scope`, `envelope`, and `loads` data.

NOTE! The `mod` keyword changes from which Module data is accessed from,
`@__MODULE__` by default.

Aggregates the properties of the allocated structures according to their `structure_type_weight` parameters.
The fenestration and ventilation properties are allocated according to `interior_air_and_furniture_weight`.
Division of internal heat gains and solar heat gains between the interior air and the different structure types
is based on the `interior_air_and_furniture_weight`, the relative surface areas of the different structure types
and their `structure_type_weights`, as well as the `internal_heat_gain_convective_fraction` parameter
for controlling the assumed convective vs radiative fractions of the internal gains.
Solar gains are divided between the interior air and structures similar to internal heat gains,
but use the `solar_heat_gain_convective_fraction` parameter for their convective vs radiative fraction instead.

Essentially, this function performs the following steps:
1. Fetch user defined thermal mass, self-discharge, temperature set-point, and heat transfer coefficient parameters.
2. Calculate the thermal mass on the `node` using [`calculate_interior_air_and_furniture_thermal_mass`](@ref) and [`calculate_structural_thermal_mass`](@ref).
3. Calculate the total heat transfer coefficient between the structures on this `node` and the interior air using [`calculate_structural_interior_heat_transfer_coefficient`](@ref).
4. Calculate the total heat transfer coefficient between the structures on this `node` and the ambient air using [`calculate_structural_exterior_heat_transfer_coefficient`](@ref).
5. Calculate the total heat transfer coefficient between the structures on this `node` and the ground using [`calculate_structural_ground_heat_transfer_coefficient`](@ref).
6. Calculate the total heat transfer coefficient through windows for this `node` using [`calculate_window_heat_transfer_coefficient`](@ref).
7. Calculate the total heat transfer coefficient of ventilation and infiltration on this `node` using [`calculate_ventilation_and_infiltration_heat_transfer_coefficient`](@ref).
8. Calculate the total heat transfer coefficient of thermal bridges for this `node` using [`calculate_linear_thermal_bridge_heat_transfer_coefficient`](@ref).
9. Fetch domestic hot water demand from `loads` for this `node`.
10. Calculate the convective internal heat gains on this `node` using [`calculate_convective_internal_heat_gains`](@ref).
11. Calculate the radiative internal heat gains on this `node` using [`calculate_radiative_internal_heat_gains`](@ref).
12. Calculate the convective solar heat gains on this `node` using [`calculate_convective_solar_gains`](@ref).
13. Calculate the radiative solar heat gains on this `node` using [`calculate_radiative_solar_gains`](@ref).
14. Return all the pieces necessary for constructing the [`BuildingNodeData`](@ref) for this `node`.

**NOTE! Linear thermal bridges are assumed to bypass any potential structural lumped-capacitance nodes,
and apply directly to the heat transfer coefficient between the interior air and the ambient air!**
**NOTE! Currently, radiative internal and solar gains are lost through the windows
in the building envelope.**
"""
function process_building_node(
    archetype::Object,
    node::Object,
    scope::ScopeData,
    envelope::EnvelopeData,
    loads::LoadsData;
    mod::Module = @__MODULE__,
)
    # Fetch the interior weight of the node.
    interior_weight = mod.interior_air_and_furniture_weight(building_node = node)

    # Fetch user defined thermal mass components of the node.
    thermal_mass_base_J_K = mod.effective_thermal_mass_base_J_K(building_node = node)
    thermal_mass_gfa_scaled_J_K =
        scope.average_gross_floor_area_m2_per_building *
        mod.effective_thermal_mass_gfa_scaling_J_m2K(building_node = node)

    # Fetch user-defined self-discharge rates.
    self_discharge_base_W_K = mod.self_discharge_rate_base_W_K(building_node = node)
    self_discharge_gfa_scaled_W_K =
        scope.average_gross_floor_area_m2_per_building *
        mod.self_discharge_rate_gfa_scaling_W_m2K(building_node = node)

    # Set temperature bounds with a potential override for indoor air temperature set points.
    if (
        interior_weight > 0 &&
        !isnothing(
            mod.indoor_air_cooling_set_point_override_K(building_archetype = archetype),
        )
    )
        maximum_temperature_K =
            mod.indoor_air_cooling_set_point_override_K(building_archetype = archetype)
    else
        maximum_temperature_K = mod.maximum_permitted_temperature_K(building_node = node)
    end
    if (
        interior_weight > 0 &&
        !isnothing(
            mod.indoor_air_heating_set_point_override_K(building_archetype = archetype),
        )
    )
        minimum_temperature_K =
            mod.indoor_air_heating_set_point_override_K(building_archetype = archetype)
    else
        minimum_temperature_K = mod.minimum_permitted_temperature_K(building_node = node)
    end

    # Fetch the user-defined heat transfer coefficients
    heat_transfer_coefficients_base_W_K = Dict(
        n => mod.heat_transfer_coefficient_base_W_K(
            building_node1 = node,
            building_node2 = n,
        ) for n in mod.building_node__building_node(building_node1 = node)
    )
    heat_transfer_coefficients_gfa_scaled_W_K = Dict(
        n =>
            scope.average_gross_floor_area_m2_per_building *
            mod.heat_transfer_coefficient_gfa_scaling_W_m2K(
                building_node1 = node,
                building_node2 = n,
            ) for n in mod.building_node__building_node(building_node1 = node)
    )

    # Calculate interior air and structural thermal masses.
    thermal_mass_interior_air_and_furniture_J_K =
        calculate_interior_air_and_furniture_thermal_mass(
            archetype,
            scope,
            interior_weight;
            mod = mod,
        )
    thermal_mass_structures_J_K =
        calculate_structural_thermal_mass(node, scope, envelope; mod = mod)

    # Calculate the heat transfer coefficients of the included structures.
    heat_transfer_coefficient_structures_interior_W_K =
        calculate_structural_interior_heat_transfer_coefficient(
            node,
            scope,
            envelope,
            interior_weight;
            mod = mod,
        )
    heat_transfer_coefficient_structures_exterior_W_K =
        calculate_structural_exterior_heat_transfer_coefficient(
            node,
            scope,
            envelope,
            interior_weight;
            mod = mod,
        )
    heat_transfer_coefficient_structures_ground_W_K =
        calculate_structural_ground_heat_transfer_coefficient(
            node,
            scope,
            envelope;
            mod = mod,
        )

    # Calculate the heat transfer coefficients for fenestration, ventilation, and thermal bridges
    heat_transfer_coefficient_windows_W_K =
        calculate_window_heat_transfer_coefficient(scope, envelope, interior_weight)
    heat_transfer_coefficient_ventilation_and_infiltration_W_K =
        calculate_ventilation_and_infiltration_heat_transfer_coefficient(
            archetype,
            scope,
            interior_weight;
            mod = mod,
        )
    heat_transfer_coefficient_thermal_bridges_W_K =
        calculate_linear_thermal_bridge_heat_transfer_coefficient(
            scope,
            envelope,
            interior_weight;
            mod = mod,
        )

    # Fetch the DHW demand, as well as internal and solar heat gains for the node.
    if mod.domestic_hot_water_demand_weight(building_node = node) > 0
        domestic_hot_water_demand_W =
            loads.domestic_hot_water_demand_W *
            mod.domestic_hot_water_demand_weight(building_node = node)
    else
        domestic_hot_water_demand_W = 0
    end

    # Calculate the total structure area for distributing the internal and solar gains.
    total_structure_area_m2 =
        sum(getfield(envelope, field).surface_area_m2 for field in fieldnames(EnvelopeData))

    # Calculate the internal heat gains for the node, first for internal air and then structures.
    internal_heat_gains_air_W = calculate_convective_internal_heat_gains(
        archetype,
        loads,
        interior_weight;
        mod = mod,
    )
    internal_heat_gains_structures_W = calculate_radiative_internal_heat_gains(
        archetype,
        node,
        envelope,
        loads,
        total_structure_area_m2;
        mod = mod,
    )

    # Calculate the solar heat gains for the node, first for internal air and then structures.
    solar_heat_gains_air_W =
        calculate_convective_solar_gains(archetype, loads, interior_weight; mod = mod)
    solar_heat_gains_structures_W = calculate_radiative_solar_gains(
        archetype,
        node,
        envelope,
        loads,
        total_structure_area_m2;
        mod = mod,
    )

    # Return all the stuff in the correct order.
    return thermal_mass_base_J_K,
    thermal_mass_gfa_scaled_J_K,
    thermal_mass_interior_air_and_furniture_J_K,
    thermal_mass_structures_J_K,
    maximum_temperature_K,
    minimum_temperature_K,
    self_discharge_base_W_K,
    self_discharge_gfa_scaled_W_K,
    heat_transfer_coefficients_base_W_K,
    heat_transfer_coefficients_gfa_scaled_W_K,
    heat_transfer_coefficient_structures_interior_W_K,
    heat_transfer_coefficient_structures_exterior_W_K,
    heat_transfer_coefficient_structures_ground_W_K,
    heat_transfer_coefficient_windows_W_K,
    heat_transfer_coefficient_ventilation_and_infiltration_W_K,
    heat_transfer_coefficient_thermal_bridges_W_K,
    domestic_hot_water_demand_W,
    internal_heat_gains_air_W,
    internal_heat_gains_structures_W,
    solar_heat_gains_air_W,
    solar_heat_gains_structures_W,
    interior_weight
end


"""
    calculate_interior_air_and_furniture_thermal_mass(
        archetype::Object,
        scope::ScopeData,
        interior_weight::Real;
        mod::Module = @__MODULE__,
    )

Calculate the effective thermal mass of interior air and furniture on `node` in [J/K].

NOTE! The `mod` keyword changes from which Module data is accessed from,
`@__MODULE__` by default.

Essentially, the calculation is based on the gross-floor area [m2] of the `archetype`
building, the assumed `effective_thermal_capacity_of_interior_air_and_furniture_J_m2K`
of the `archetype`, and the assumed `interior_air_and_furniture_weight` of the `node`.
```math
C_\\text{int,n} = w_\\text{int,n} c_\\text{int,gfa} A_\\text{gfa}
```
"""
function calculate_interior_air_and_furniture_thermal_mass(
    archetype::Object,
    scope::ScopeData,
    interior_weight::Real;
    mod::Module = @__MODULE__,
)
    scope.average_gross_floor_area_m2_per_building *
    mod.effective_thermal_capacity_of_interior_air_and_furniture_J_m2K(
        building_archetype = archetype,
    ) *
    interior_weight
end


"""
    calculate_structural_thermal_mass(
        node::Object,
        scope::ScopeData,
        envelope::EnvelopeData;
        mod::Module = @__MODULE__,
    )

Calculate the effective thermal mass of the structures on `node` in [J/K].

NOTE! The `mod` keyword changes from which Module data is accessed from,
`@__MODULE__` by default.

Essentially, sums the total effective thermal mass of the structures
attributed to this node, accounting for their `structure_type_weight`s.
```math
C_\\text{n,str} = \\sum_{\\text{st} \\in n} w_\\text{n,st} c_\\text{st} A_\\text{st}
```
"""
function calculate_structural_thermal_mass(
    node::Object,
    scope::ScopeData,
    envelope::EnvelopeData;
    mod::Module = @__MODULE__,
)
    reduce(
        +,
        scope.structure_data[st].effective_thermal_mass_J_m2K *
        getfield(envelope, st.name).surface_area_m2 *
        mod.structure_type_weight(building_node = node, structure_type = st) for
        st in mod.building_node__structure_type(building_node = node);
        init = 0,
    )
end


"""
    calculate_structural_interior_heat_transfer_coefficient(
        node::Object,
        scope::ScopeData,
        envelope::EnvelopeData,
        interior_weight::Real;
        mod::Module = @__MODULE__,
    )

Calculate the total interior heat transfer coefficient of the structures in `node` in [W/K].

NOTE! The `mod` keyword changes from which Module data is accessed from,
`@__MODULE__` by default.

Essentially, initializes the total heat transfer coefficient between
the interior air node and the node containing the structures.
```math
H_\\text{int,n} = \\begin{cases}
(1 - w_\\text{int,n}) \\sum_{\\text{st} \\in n} w_\\text{n,st} U_\\text{int,st} A_\\text{st} \\qquad \\text{st} \\notin \\text{internal structures} \\\\
(1 - w_\\text{int,n}) \\sum_{\\text{st} \\in n} w_\\text{n,st} ( U_\\text{int,st} + U_\\text{ext,st} ) A_\\text{st} \\qquad \\text{st} \\notin \\text{internal structures}
\\end{cases}
```

**NOTE!** For `is_internal` structures, the U-values for both the interior and
*exterior* surface are used, effectively accounting for both sides of interior
structures.
**NOTE!** If an external structure is lumped together with the
interior air node, it's internal heat transfer coefficient is attributed to its
exterior heat transfer instead, *but using this feature is not recommended*!
"""
function calculate_structural_interior_heat_transfer_coefficient(
    node::Object,
    scope::ScopeData,
    envelope::EnvelopeData,
    interior_weight::Real;
    mod::Module = @__MODULE__,
)
    reduce(
        +,
        (
            scope.structure_data[st].internal_U_value_to_structure_W_m2K + (
                mod.is_internal(structure_type = st) ?
                scope.structure_data[st].external_U_value_to_ambient_air_W_m2K : 0
            )
        ) *
        (1 - interior_weight) * # If an external structure is lumped together with air node, it's internal heat coefficient adds to the external heat transfer.
        getfield(envelope, st.name).surface_area_m2 *
        mod.structure_type_weight(building_node = node, structure_type = st) for
        st in mod.building_node__structure_type(building_node = node);
        init = 0,
    )
end


"""
    calculate_structural_exterior_heat_transfer_coefficient(
        node::Object,
        scope::ScopeData,
        envelope::EnvelopeData,
        interior_weight::Real
    )

Calculate the total exterior heat transfer coefficient of the structures in `node` in [W/K].

NOTE! The `mod` keyword changes from which Module data is accessed from,
`@__MODULE__` by default.

Essentially, calculates the total heat transfer coefficient between
the ambient air and the node containing the structures.
Internal structures have no exterior heat transfer coefficient,
as they are assumed to not be part of the building envelope.
```math
H_\\text{ext,n} = \\sum_{\\text{st} \\in n} w_\\text{st} \\left( \\frac{1}{U_\\text{ext,st}} + \\frac{w_{int,n}}{U_\\text{int,st}} \\right)^{-1} A_\\text{st}
```

**NOTE!** If an external structure is lumped together with the
interior air node, it's internal heat transfer coefficient is attributed to its
exterior heat transfer as well, *but using this feature is not recommended*!
"""
function calculate_structural_exterior_heat_transfer_coefficient(
    node::Object,
    scope::ScopeData,
    envelope::EnvelopeData,
    interior_weight::Real;
    mod::Module = @__MODULE__,
)
    reduce(
        +,
        1 / (
            1 / scope.structure_data[st].external_U_value_to_ambient_air_W_m2K +
            interior_weight / # If an external structure is lumped together with air, it's internal heat coefficient adds to the external heat transfer.
            scope.structure_data[st].internal_U_value_to_structure_W_m2K
        ) *
        getfield(envelope, st.name).surface_area_m2 *
        mod.structure_type_weight(building_node = node, structure_type = st) for
        st in mod.building_node__structure_type(building_node = node) if
        !(mod.is_internal(structure_type = st));
        init = 0,
    )
end


"""
    calculate_structural_ground_heat_transfer_coefficient(
        node::Object,
        scope::ScopeData,
        envelope::EnvelopeData;
        mod::Module = @__MODULE__,
    )

Calculate the total ground heat transfer coefficient of the structures in `node` in [W/K].

NOTE! The `mod` keyword changes from which Module data is accessed from,
`@__MODULE__` by default.

Essentially, calculates the total heat transfer coefficient between
the ground and the node containing the structures.
```math
H_\\text{grn,n} = \\sum_{\\text{st} \\in n} w_\\text{st} U_\\text{grn,st} A_\\text{st}
```
"""
function calculate_structural_ground_heat_transfer_coefficient(
    node::Object,
    scope::ScopeData,
    envelope::EnvelopeData;
    mod::Module = @__MODULE__,
)
    reduce(
        +,
        scope.structure_data[st].external_U_value_to_ground_W_m2K *
        getfield(envelope, st.name).surface_area_m2 *
        mod.structure_type_weight(building_node = node, structure_type = st) for
        st in mod.building_node__structure_type(building_node = node);
        init = 0,
    )
end


"""
    calculate_window_heat_transfer_coefficient(
        scope::ScopeData,
        envelope::EnvelopeData,
        interior_weight::Real
    )
Calculate the total heat transfer coefficient for windows.

Windows are assumed to transfer heat directly between the indoor air and the ambient air.
```math
H_\\text{w,n} = w_\\text{int,n} U_\\text{w} A_\\text{w}
```
"""
function calculate_window_heat_transfer_coefficient(
    scope::ScopeData,
    envelope::EnvelopeData,
    interior_weight::Real,
)
    scope.window_U_value_W_m2K * envelope.window.surface_area_m2 * interior_weight
end


"""
    calculate_ventilation_and_infiltration_heat_transfer_coefficient(
        archetype::Object,
        scope::ScopeData,
        interior_weight::Real;
        mod::Module = @__MODULE__,
    )
    
Calculate ventilation and infiltration heat transfer coefficient.

NOTE! The `mod` keyword changes from which Module data is accessed from,
`@__MODULE__` by default.

Ventilation and infiltration are assumed to transfer heat directly between the
interior and ambient air.
Loosely based on EN ISO 52016-1:2017 6.5.10.1.
```math
H_\\text{ven,n} = w_\\text{int,n} A_\\text{gfa} h_\\text{room} \\rho_\\text{air} \\frac{ (1 - \\eta_\\text{hru}) r_\\text{ven} + r_\\text{inf}}{3600}
```
"""
function calculate_ventilation_and_infiltration_heat_transfer_coefficient(
    archetype::Object,
    scope::ScopeData,
    interior_weight::Real;
    mod::Module = @__MODULE__,
)
    mod.room_height_m(building_archetype = archetype) *
    scope.average_gross_floor_area_m2_per_building *
    mod.volumetric_heat_capacity_of_interior_air_J_m3K(building_archetype = archetype) /
    3600 *
    (
        scope.ventilation_rate_1_h * (1 - scope.HRU_efficiency) +
        scope.infiltration_rate_1_h
    ) *
    interior_weight
end


"""
    calculate_linear_thermal_bridge_heat_transfer_coefficient(
        scope::ScopeData,
        envelope::EnvelopeData,
        interior_weight::Real;
        mod::Module = @__MODULE__,
    )

Calculate linear thermal bridge heat transfer coefficient.

NOTE! The `mod` keyword changes from which Module data is accessed from,
`@__MODULE__` by default.

Linear thermal bridges are assumed to bypass the temperature node within the structure,
and act as direct heat transfer between the indoor air and ambient conditions.
```math
H_{\\Psi,n} = w_\\text{int,n} \\sum_\\text{st} l_\\text{st} \\Psi_\\text{st}
```
"""
function calculate_linear_thermal_bridge_heat_transfer_coefficient(
    scope::ScopeData,
    envelope::EnvelopeData,
    interior_weight::Real;
    mod::Module = @__MODULE__,
)
    reduce(
        +,
        scope.structure_data[st].linear_thermal_bridges_W_mK *
        getfield(envelope, st.name).linear_thermal_bridge_length_m for
        st in mod.structure_type();
        init = 0,
    ) * interior_weight
end


"""
    calculate_convective_internal_heat_gains(
        archetype::Object
        loads::LoadsData,
        interior_weight::Real;
        mod::Module = @__MODULE__,
    )

Calculate the convective internal heat gains on the `node` in [W].

NOTE! The `mod` keyword changes from which Module data is accessed from,
`@__MODULE__` by default.

Essentially, takes the given internal heat gain profile in `loads` and
multiplies it with the share of interior air on this `node` as well as the
assumed convective fraction of internal heat gains.
```math
\\Phi_\\text{int,conv,n} = w_\\text{int,n} f_\\text{int,conv} \\Phi_\\text{int}
```
"""
function calculate_convective_internal_heat_gains(
    archetype::Object,
    loads::LoadsData,
    interior_weight::Real;
    mod::Module = @__MODULE__,
)
    loads.internal_heat_gains_W *
    interior_weight *
    mod.internal_heat_gain_convective_fraction(building_archetype = archetype)
end


"""
    calculate_radiative_internal_heat_gains(
        archetype::Object,
        node::Object,
        envelope::EnvelopeData,
        loads::LoadsData,
        total_structure_area_m2::Real;
        mod::Module = @__MODULE__,
    )

Calculate the radiative internal heat gains on the `node` in [W].

NOTE! The `mod` keyword changes from which Module data is accessed from,
`@__MODULE__` by default.

Essentially, takes the given internal heat gain profile in `loads`
and multiplies it with the assumed radiative fraction of internal heat gains.
The radiative heat gains are assumed to be distributed across the structures 
simply based on their relative surface areas.
```math
\\Phi_\\text{int,rad,n} = (1 - f_\\text{int,conv}) \\frac{\\sum_{\\text{st} \\in n} A_\\text{st}}{\\sum_{\\text{st}} A_\\text{st}} \\Phi_\\text{int}
```

**NOTE!** Currently, radiative internal heat gains are lost through windows!
"""
function calculate_radiative_internal_heat_gains(
    archetype::Object,
    node::Object,
    envelope::EnvelopeData,
    loads::LoadsData,
    total_structure_area_m2::Real;
    mod::Module = @__MODULE__,
)
    loads.internal_heat_gains_W *
    (1 - mod.internal_heat_gain_convective_fraction(building_archetype = archetype)) *
    reduce(
        +,
        getfield(envelope, st.name).surface_area_m2 / total_structure_area_m2 for
        st in mod.building_node__structure_type(building_node = node);
        init = 0,
    )
end


"""
    calculate_convective_solar_gains(
        archetype::Object
        loads::LoadsData,
        interior_weight::Real;
        mod::Module = @__MODULE__,
    )

Calculate the convective solar heat gains on the `node` in [W].

NOTE! The `mod` keyword changes from which Module data is accessed from,
`@__MODULE__` by default.

Essentially, takes the given solar heat gain profile in `loads` and
multiplies it with the share of interior air on this `node` as well as the
assumed convective fraction of solar heat gains.
```math
\\Phi_\\text{sol,conv,n} = w_\\text{int,n} f_\\text{sol,conv} \\Phi_\\text{sol}
```
"""
function calculate_convective_solar_gains(
    archetype::Object,
    loads::LoadsData,
    interior_weight::Real;
    mod::Module = @__MODULE__,
)
    loads.solar_heat_gains_W *
    interior_weight *
    mod.solar_heat_gain_convective_fraction(building_archetype = archetype)
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

Calculate the radiative solar heat gains on the `node` in [W].

NOTE! The `mod` keyword changes from which Module data is accessed from,
`@__MODULE__` by default.

Essentially, takes the given solar heat gain profile in `loads`
and multiplies it with the assumed radiative fraction of solar heat gains.
The radiative heat gains are assumed to be distributed across the structures 
simply based on their relative surface areas.
```math
\\Phi_\\text{sol,rad,n} = (1 - f_\\text{sol,conv}) \\frac{\\sum_{\\text{st} \\in n} A_\\text{st}}{\\sum_{\\text{st}} A_\\text{st}} \\Phi_\\text{sol}
```

**NOTE!** Currently, radiative solar heat gains are lost through windows!
"""
function calculate_radiative_solar_gains(
    archetype::Object,
    node::Object,
    envelope::EnvelopeData,
    loads::LoadsData,
    total_structure_area_m2::Real;
    mod::Module = @__MODULE__,
)
    loads.solar_heat_gains_W *
    (1 - mod.solar_heat_gain_convective_fraction(building_archetype = archetype)) *
    reduce(
        +,
        getfield(envelope, st.name).surface_area_m2 / total_structure_area_m2 for
        st in mod.building_node__structure_type(building_node = node);
        init = 0,
    )
end


"""
    create_abstract_node_network(
        building_node_network::BuildingNodeNetwork,
        weather::WeatherData
    )

Process a `BuildingNodeNetwork` into an `AbstractNodeNetwork`.

The `AbstractNodeNetwork` is a useful step for creating model-agnostic input
for multiple large-scale energy system models.
`weather` is required to account for ambient temperatures.
"""
function create_abstract_node_network(
    building_node_network::BuildingNodeNetwork,
    weather::WeatherData,
)
    Dict(
        node => AbstractNode(building_node_network, node, weather) for
        node in keys(building_node_network)
    )::AbstractNodeNetwork
end


"""
    process_abstract_node(
        building_node_network::BuildingNodeNetwork,
        node::Object,
        weather::WeatherData
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
5. Sum together the internal heat gains, solar gains, DHW demand, as well as the impact of ambient temperatures.
6. Return the components required for constructing an [`AbstractNode`](@ref).

**NOTE!** The ambient temperatures are accounted for via a combination of `self_discharge_coefficient_W_K`
and `external_load`, instead of  `heat_transfer_coefficients_W_K` on any ambient temperature nodes.
This is because not all large-scale energy system models support ambient temperatures as separate parameters,
whereas self-discharge and external loads are almost always supported.
The principle is illustrated by the equation below:
```math
\\Phi_\\text{ambient heat losses} = H_\\text{ext}(T_\\text{ambient} - T_\\text{internal}) \\\\
= H_\\text{ext} T_\\text{ambient} - H_\\text{ext} T_\\text{internal} \\\\
= \\Phi_\\text{ambient} - \\Phi_\\text{self discharge}
```

**NOTE!** All heat transfer coefficients are assumed to be symmetrical!
**NOTE!** All [`AbstractNode`](@ref)s are given `1e-9 Wh/K` thermal mass to avoid
singularities when solving the temperature dynamics and heat demand later on.
"""
function process_abstract_node(
    building_node_network::BuildingNodeNetwork,
    node::Object,
    weather::WeatherData,
)
    # Convenience access to the `BuildingNodeData`.
    node_data = building_node_network[node]

    # Total thermal mass of the node, in Wh/K for better scaling in energy system models.
    thermal_mass_Wh_K =
        (
            node_data.thermal_mass_base_J_K +
            node_data.thermal_mass_gfa_scaled_J_K +
            node_data.thermal_mass_interior_air_and_furniture_J_K +
            node_data.thermal_mass_structures_J_K
        ) / 3600 + 1e-9 # Token thermal mass always required to avoid singularities in the dynamics matrix

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
    # First, connection to interior air.
    heat_transfer_coefficients_W_K = Dict(
        n =>
            node_data.heat_transfer_coefficient_structures_interior_W_K *
            building_node_network[n].interior_air_and_furniture_weight for
        n in keys(building_node_network)
    )
    # Force symmetry.
    mergewith!(
        +,
        heat_transfer_coefficients_W_K,
        Dict(
            n =>
                building_node_network[n].heat_transfer_coefficient_structures_interior_W_K *
                node_data.interior_air_and_furniture_weight for
            n in keys(building_node_network)
        ),
    )
    # Then updated with user-defined heat-transfer coefficients
    mergewith!(
        +,
        heat_transfer_coefficients_W_K,
        node_data.heat_transfer_coefficients_base_W_K,
    )
    mergewith!(
        +,
        heat_transfer_coefficients_W_K,
        node_data.heat_transfer_coefficients_gfa_scaled_W_K,
    )
    # Force symmetry.
    for n in keys(building_node_network)
        heat_transfer_coefficients_W_K[n] += (
            get(building_node_network[n].heat_transfer_coefficients_base_W_K, node, 0) +
            get(
                building_node_network[n].heat_transfer_coefficients_gfa_scaled_W_K,
                node,
                0,
            )
        )
    end
    # And filter out zero heat transfer coefficients.
    filter!(pair -> pair[2] != 0, heat_transfer_coefficients_W_K)

    # External load accounting for heat transfer with ambient conditions.
    external_load =
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
        node_data.solar_heat_gains_air_W +
        node_data.solar_heat_gains_structures_W - node_data.domestic_hot_water_demand_W

    # Return the properties of interest in the correct order for `AbstractNode`.
    return thermal_mass_Wh_K,
    self_discharge_coefficient_W_K,
    heat_transfer_coefficients_W_K,
    external_load,
    node_data.minimum_temperature_K,
    node_data.maximum_temperature_K
end
