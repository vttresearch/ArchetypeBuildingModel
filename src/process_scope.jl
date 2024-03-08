#=
    process_scope.jl

Contains functions for processing the data for `building_scope` database objects.
=#


"""
    process_building_stock_scope(scope::Object; mod::Module = @__MODULE__)

Process the aggregated `building_stock_statistic` properties of a `scope`.

NOTE! The `mod` keyword changes from which Module data is accessed from,
`@__MODULE__` by default.

Returns the gross-floor area weighted total `number_of_buildings` and
`average_gross_floor_area_m2_per_building` defined by the `building_scope` and its associated relationships.
Furthermore, returns the gross-floor area weights for the `building_type__building_period__location_id` useful
for weighting the structural and fenestration/ventilation data.

Essentially, performs the following steps:
1. Filter out irrelevant building stock statistics and check there are any left.
2. Calculate the gross-floor area weights using [`calculate_gross_floor_area_weights`](@ref).
3. Aggregate the gross-floor area weights using [`aggregate_gfa_weights`](@ref).
4. Return the necessary pieces to construct a [`ScopeData`](@ref).
"""
function process_building_stock_scope(scope::Object; mod::Module=@__MODULE__)
    # Find the relevant building stock statistics with all `building_periods`.
    relevant_building_stock_statistics = mod.building_stock_statistics(
        building_stock=mod.building_scope__building_stock(building_scope=scope),
        building_type=mod.building_scope__building_type(building_scope=scope),
        location_id=mod.building_scope__location_id(building_scope=scope),
        heat_source=mod.building_scope__heat_source(building_scope=scope);
        _compact=false
    )

    # Find the relevant `building_period` weights as limited by `scope_period_start_year` and `scope_period_end_year`
    building_period_weights = Dict(
        bp => _building_period_weight(bp, scope; mod=mod) for
        bp in unique(getfield.(relevant_building_stock_statistics, :building_period))
    )
    filter!(pair -> pair[2] > 0, building_period_weights)

    # Filter the relevant building stock statistics to only contain the relevant `building_periods`.
    relevant_building_stock_statistics = filter( # We can't use `filter!` here without affecting the underlying `building_stock_statistics`
        namtup -> namtup.building_period in keys(building_period_weights),
        relevant_building_stock_statistics,
    )

    # Check that some statistics are found, error otherwise.
    if isempty(relevant_building_stock_statistics)
        @error """
        No relevant building stock statistics found for `building_scope` `$(scope)`!
        Please check the definitions, e.g. `scope_period_end_year` and `scope_period_start_year` parameters.
        """
    end

    # Calculate the gross-floor area weights and other useful indicators.
    gross_floor_area_weights,
    total_weighted_number_of_buildings,
    total_weighted_gross_floor_area = calculate_gross_floor_area_weights(
        scope,
        relevant_building_stock_statistics,
        building_period_weights;
        mod=mod
    )

    # Calculate the aggregated weights for ventilation/infiltration and weather data processing.
    aggregated_gfa_weights, location_id_gfa_weights =
        aggregate_gfa_weights(gross_floor_area_weights)

    return total_weighted_gross_floor_area / total_weighted_number_of_buildings,
    total_weighted_number_of_buildings,
    gross_floor_area_weights,
    aggregated_gfa_weights,
    location_id_gfa_weights
end


"""
    _building_period_weight(building_period::Object, building_scope::Object; mod::Module = @__MODULE__)

Calculate the weight of a `building_period` within a `building_scope`.

NOTE! The `mod` keyword changes from which Module data is accessed from,
`@__MODULE__` by default.

Essentially, represents whether the [building\\_period](@ref) `bp` is contained
within [scope\\_period\\_start\\_year](@ref) -- [scope\\_period\\_end\\_year](@ref)
in its entirety [1], only partially (0,1), or if at all [0].
```math
w_\\text{bp} = \\text{max}\\left( \\text{min} \\left( \\frac{\\text{end}_\\text{scope} - \\text{start}_\\text{bp}}{\\text{end}_\\text{bp} - \\text{start}_\\text{bp}}, \\frac{\\text{end}_\\text{bp} - \\text{start}_\\text{scope}}{\\text{end}_\\text{bp} - \\text{start}_\\text{bp}}, 1 \\right), 0 \\right)
```
"""
function _building_period_weight(period::Object, scope::Object; mod::Module=@__MODULE__)
    # Fetch period and scope years.
    bp_start = mod.period_start(building_period=period)
    bp_end = mod.period_end(building_period=period)
    scope_start = mod.scope_period_start_year(building_scope=scope)
    scope_end = mod.scope_period_end_year(building_scope=scope)

    # Check that the years are sensible.
    if bp_end < bp_start
        @error """
        End year cannot be less than start year for `building_period` `$(period)`!
        """
    elseif scope_end < scope_start
        @error """
        End year cannot be less than start year for `building_scope` `$(scope)`!
        """
    end

    # Calculate the weight.
    weight = max(
        min(
            (scope_end - bp_start) / (bp_end - bp_start),
            (bp_end - scope_start) / (bp_end - bp_start),
            1,
        ),
        0,
    )
    isnan(weight) ? 1 : weight
end


"""
    calculate_gross_floor_area_weights(
        scope::Object,
        relevant_building_stock_statistics::Vector,
        building_period_weights::Dict{Object,T} where T <: Real;
        mod::Module = @__MODULE__,
    )

Calculate the gross-floor area weights for the `relevant_building_stock_statistics`.

NOTE! The `mod` keyword changes from which Module data is accessed from,
`@__MODULE__` by default.

Essentially, returns a normalized weight for each relevant entry in the input
data statistics, corresponding to how impactful it is when aggregating
the properties for the given `scope`.
The weights of the individual dimensions are assumed to be one,
unless otherwise specified in [The `building_scope` definition](@ref).
Also returns the total number of buildings and
the total weighted gross-floor area included in the `building_scope`.
```math
w_\\text{bs,bt,bp,lid,hs} = \\frac{w_\\text{bs} w_\\text{bt} w_\\text{bp} w_\\text{lid} w_\\text{hs} n_\\text{bs,bt,bp,lid,hs} A_\\text{gfa,bs,bt,bp,lid,hs}}{\\sum_{bs,bt,bp,lid,hs} w_\\text{bs,bt,bp,lid,hs}}
```
where `w_bs` is the [building\\_stock\\_weight](@ref),
`w_bt` is the [building\\_type\\_weight](@ref),
`w_bp` is the [`_building_period_weight`](@ref),
`w_lid` is the [location\\_id\\_weight](@ref),
`w_hs` is the [heat\\_source\\_weight](@ref),
`n_bs,bt,bp,lid,hs` is the [number\\_of\\_buildings](@ref),
and `A_gfa` is the [average\\_gross\\_floor\\_area\\_m2\\_per\\_building](@ref).
"""
function calculate_gross_floor_area_weights(
    scope::Object,
    relevant_building_stock_statistics::Vector,
    building_period_weights::Dict{Object,T} where {T<:Real};
    mod::Module=@__MODULE__
)
    # Initialize gross-floor area weights by calculating the weighted number of buildings.
    gross_floor_area_weights = [
        (bs, bt, bp, lid, hs) =>
            mod.number_of_buildings(
                building_stock=bs,
                building_type=bt,
                building_period=bp,
                location_id=lid,
                heat_source=hs,
            ) *
            mod.building_stock_weight(building_scope=scope, building_stock=bs) *
            mod.building_type_weight(building_scope=scope, building_type=bt) *
            building_period_weights[bp] *
            mod.location_id_weight(building_scope=scope, location_id=lid) *
            mod.heat_source_weight(building_scope=scope, heat_source=hs) for
        (bs, bt, bp, lid, hs) in relevant_building_stock_statistics
    ]
    # Record total weighted number of buildings and check that its non-zero.
    total_weighted_number_of_buildings = sum(getfield.(gross_floor_area_weights, 2))
    if isapprox(total_weighted_number_of_buildings, 0)
        error(
            """
            `total_weighted_number_of_buildings` of `building_scope=$(scope)` is approximately zero!
            Check that there's meaningful data contained by `$(scope)`.
            """,
        )
    end

    # Next, multiply the weighted number of buildings by the average gross-floor area to obtain total GFA.
    gross_floor_area_weights = [
        (bs, bt, bp, lid, hs) =>
            val * mod.average_gross_floor_area_m2_per_building(
                building_stock=bs,
                building_type=bt,
                building_period=bp,
                location_id=lid,
                heat_source=hs,
            ) for ((bs, bt, bp, lid, hs), val) in gross_floor_area_weights
    ]
    # Calculate total weighted GFA, and use it to normalize the weights.
    total_weighted_gross_floor_area = sum(getfield.(gross_floor_area_weights, 2))
    gross_floor_area_weights = Dict(
        (bs, bt, bp, lid, hs) => val / total_weighted_gross_floor_area for
        ((bs, bt, bp, lid, hs), val) in gross_floor_area_weights
    )

    # Filter out useless entries and check that the total weight is approximately 1.
    filter!(pair -> pair[2] > 0, gross_floor_area_weights)
    if !isapprox(sum(values(gross_floor_area_weights)), 1)
        error("""
              Gross-floor area weights of `$(scope)` don't sum to 1!
              Got $(sum(values(gross_floor_area_weights))) instead.
              """)
    end

    # Return the final gross floor area weights and weighted total number of buildings.
    return gross_floor_area_weights,
    total_weighted_number_of_buildings,
    total_weighted_gross_floor_area
end


"""
    aggregate_gfa_weights(gross_floor_area_weights::Dict{NTuple{5,Object},Float64})

Aggregate `gross_floor_area_weights` for ventilation/infiltration and weather data processing.

Essentially, sums over the gross-floor area weights of unused dimensions to
produce reduced sets of weights.
```math
w_\\text{bt,bp,lid} = \\sum_\\text{bp,hs} w_\\text{bs,bt,bp,lid,hs} \\\\
w_\\text{lid} = \\sum_\\text{bt,bp} w_\\text{bt,bp,lid},
```
where `w_bs,bt,bp,lid,hs` are the full gross-floor area weights calculated
using [`calculate_gross_floor_area_weights`](@ref).
"""
function aggregate_gfa_weights(gross_floor_area_weights::Dict{NTuple{5,Object},Float64})
    # Aggregate the GFA-weights to match the dimensions of the ventilation and fenestration (and structural) data.
    aggregated_gfa_weights = Dict{NTuple{3,Object},Float64}(
        (bt, bp, lid) => 0.0
        for (bs, bt, bp, lid, hs) in keys(gross_floor_area_weights)
    )
    for ((bs, bt, bp, lid, hs), weight) in gross_floor_area_weights
        aggregated_gfa_weights[(bt, bp, lid)] += weight
    end
    if !isapprox(sum(values(aggregated_gfa_weights)), 1)
        error("""
              Aggregated gross-floor area weights of `$(scope)` don't sum to 1!
              Got $(sum(values(aggregated_gfa_weights))) instead.
              """)
    end

    # Aggregate the GFA-weights for potential weather data processing.
    location_id_gfa_weights = Dict{Object,Float64}(
        lid => 0.0
        for (bt, bp, lid) in keys(aggregated_gfa_weights)
    )
    for ((bt, bp, lid), weight) in aggregated_gfa_weights
        location_id_gfa_weights[lid] += weight
    end
    if !isapprox(sum(values(location_id_gfa_weights)), 1)
        error(
            """
            `location_id` aggregated gross-floor area weights of `$(scope)` don't sum to 1!
            Got $(sum(values(location_id_gfa_weights))) instead.
            """,
        )
    end

    # Filter and return the weights
    return filter!(pair -> pair[2] != 0, aggregated_gfa_weights),
    filter!(pair -> pair[2] != 0, location_id_gfa_weights)
end


"""
    process_ventilation_and_fenestration_scope(
        aggregated_gfa_weights::Dict{NTuple{3,Object},Float64};
        mod::Module = @__MODULE__,
    )

Aggregate `ventilation_and_fenestration_statistics` using gross-floor area weights.

NOTE! The `mod` keyword changes from which Module data is accessed from,
`@__MODULE__` by default.

Essentially, this function calculates the aggregate properties of
ventilation and fenestration based on the aggregated gross-floor area weights
produced by the [`aggregate_gfa_weights`](@ref) function.
Returns the gross-floor area weighted average HRU-efficiency, infiltration rate [1/h],
ventilation rate [1/h], window U-values [W/m2K],
and window total normal solar energy transmittance.
```math
p_\\text{avg} = \\sum_{\\text{bt,bp,lid}} w_\\text{bt,bp,lid} p_\\text{bt,bp,lid}
```
where `p` is any of the above listed properties to be aggregated. 
"""
function process_ventilation_and_fenestration_scope(
    aggregated_gfa_weights::Dict{NTuple{3,Object},Float64};
    mod::Module=@__MODULE__
)
    # Calculate the aggregated ventilation and fenestration parameters.
    hru_efficiency, inf_rate, sol_transm, ven_rate, win_U = [
        sum(
            param(building_type=bt, building_period=bp, location_id=lid) * weight
            for ((bt, bp, lid), weight) in aggregated_gfa_weights
        ) for param in [
            mod.HRU_efficiency,
            mod.infiltration_rate_1_h,
            mod.total_normal_solar_energy_transmittance,
            mod.ventilation_rate_1_h,
            mod.window_U_value_W_m2K,
        ]
    ]
    return hru_efficiency, inf_rate, sol_transm, ven_rate, win_U, aggregated_gfa_weights
end


"""
    process_structure_scope(
        aggregated_gfa_weights::Dict{NTuple{3,Object},Float64};
        mod::Module = @__MODULE__,
    )

Aggregate `structure_statistics` using gross-floor area weights.

NOTE! The `mod` keyword changes from which Module data is accessed from,
`@__MODULE__` by default.

Essentially, this function calculates the aggregate properties of structures
based on the aggregated gross floor area weights produced by
the [`aggregate_gfa_weights`](@ref) function.
Returns the gross-floor area weighted average structural properties
as [`StructureData`](@ref) for different `structure_type`s.
These include the design U-values [W/m2K], effective thermal mass [J/m2K],
linear thermal bridges [W/mK], and and U-values for heat transfer between
the building interior and the structure, the structure and ambient air,
the structure and ground, as well as the total U-value through the structure.
```math
p_\\text{avg} = \\sum_{\\text{bt,bp,lid}} w_\\text{bt,bp,lid} p_\\text{bt,bp,lid}
```
where `p` is any of the above listed properties to be aggregated.
"""
function process_structure_scope(
    aggregated_gfa_weights::Dict{NTuple{3,Object},Float64};
    mod::Module=@__MODULE__
)
    return Dict{Object,StructureData}(
        st => StructureData(
            st,
            [
                sum(
                    param(
                        building_type=bt,
                        building_period=bp,
                        location_id=lid,
                        structure_type=st,
                    ) * weight for ((bt, bp, lid), weight) in aggregated_gfa_weights
                ) for param in [
                    mod.design_U_value_W_m2K,
                    mod.effective_thermal_mass_J_m2K,
                    mod.external_U_value_to_ambient_air_W_m2K,
                    mod.external_U_value_to_ground_W_m2K,
                    mod.internal_U_value_to_structure_W_m2K,
                    mod.linear_thermal_bridges_W_mK,
                    mod.total_U_value_W_m2K,
                ]
            ]...
        )
        for st in mod.structure_type()
    )
end
