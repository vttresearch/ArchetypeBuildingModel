#=
    process_loads.jl

Contains functions for processing the domestic hot water demand, internal heat gains,
and solar gains for the archetype buildings.
=#


"""
    process_building_loads(
        archetype::Object,
        scope::ScopeData;
        mod::Module = @__MODULE__,
    )

Calculate the domestic hot water demand and internal heat gains for the archetype building.

NOTE! The `mod` keyword changes from which Module data is accessed from,
`@__MODULE__` by default.

Essentially, performs the following steps:
1. Finds the `building_loads` object corresponding to the `building_archetype`.
2. Calculates total domestic hot water (DHW) demand using [`calculate_total_dhw_demand`](@ref).
3. Calculates total internal heat loads using [`calculate_total_internal_heat_loads`](@ref).
4. Returns the calculated DHW demand and internal heat gains.
"""
function process_building_loads(
    archetype::Object,
    scope::ScopeData;
    mod::Module=@__MODULE__
)
    # Find the `building_loads` connected to the `archetype`
    loads = only(mod.building_archetype__building_loads(building_archetype=archetype))

    # Calculate loads
    dhw_demand = calculate_total_dhw_demand(loads, scope; mod=mod)
    internal_gains = calculate_total_internal_heat_loads(loads, scope; mod=mod)
    return dhw_demand,
    internal_gains
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
    mod::Module=@__MODULE__
)
    mod.domestic_hot_water_demand_base_W(building_loads=loads) +
    scope.average_gross_floor_area_m2_per_building *
    mod.domestic_hot_water_demand_gfa_scaling_W_m2(building_loads=loads)
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
    mod::Module=@__MODULE__
)
    mod.internal_heat_loads_base_W(building_loads=loads) +
    scope.average_gross_floor_area_m2_per_building *
    mod.internal_heat_loads_gfa_scaling_W_m2(building_loads=loads)
end
