#=
    process_abstract_node.jl

Contains functions for processing the properties of
abstract lumped-capacitance thermal nodes, in preparation for conversion to
large-scale energy system model input.
=#

"""
    create_abstract_node_network(
        building_node_network::BuildingNodeNetwork,
    )

Process a `BuildingNodeNetwork` into an `AbstractNodeNetwork`.

TODO: Revise documentation, rename AbstractNode to FlexibilityNode?

The `AbstractNodeNetwork` is a useful step for creating model-agnostic input
for multiple large-scale energy system models.
"""
function create_abstract_node_network(
    building_node_network::BuildingNodeNetwork,
)
    Dict(
        node => AbstractNode(building_node_network, node) for
        node in keys(building_node_network)
    )::AbstractNodeNetwork
end


"""
    process_abstract_node(
        building_node_network::BuildingNodeNetwork,
        node::Object
    )

Calculate the properties of an [`AbstractNode`](@ref) corresponding to the `node` in the `building_node_network`.

TODO: Revise documentation, rename AbstractNode to FlexibilityNode?

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

**NOTE!** All heat transfer coefficients are assumed to be symmetrical!
**NOTE!** All [`AbstractNode`](@ref)s are given `1e-9 Wh/K` thermal mass to avoid
singularities when solving the temperature dynamics and heat demand later on.
"""
function process_abstract_node(
    building_node_network::BuildingNodeNetwork,
    node::Object,
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
        ) / 3600 # TODO: is this necessary? + 1e-9 # Token thermal mass always required to avoid singularities in the dynamics matrix

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
    #= TODO: Temove as obsolete?
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
        node_data.solar_heat_gains_air_W +
        node_data.solar_heat_gains_structures_W +
        node_data.solar_heat_gains_envelope_W - node_data.radiative_envelope_sky_losses_W -
        node_data.domestic_hot_water_demand_W
    =#

    # Return the properties of interest in the correct order for `AbstractNode`.
    return thermal_mass_Wh_K,
    self_discharge_coefficient_W_K,
    heat_transfer_coefficients_W_K,
    node_data.minimum_temperature_K,
    node_data.maximum_temperature_K
end
