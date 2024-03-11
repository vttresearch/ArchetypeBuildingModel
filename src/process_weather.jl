#=
    process_weather.jl

Contains functions for processing weather data.
=#


"""
    process_weather(
        archetype::Object,
        scope_data::ScopeData,
        envelope_data::EnvelopeData,
        building_nodes::BuildingNodeNetwork;
        ignore_year::Bool=false,
        repeat::Bool=false,
        save_layouts::Bool=true,
        resampling::Int=5,
        mod::Module=@__MODULE__,
        realization::Symbol=:realization
    )

Process weather data for the [`WeatherData`](@ref) constructor.

NOTE! The `mod` keyword changes from which Module data is accessed from,
`@__MODULE__` by default. The `realization` scenario is required for processing
of effective ground temperature.

Essentially, performs the following steps:
1. Call [`create_building_weather`](@ref) to handle weather and demand processing using `ArBuWe.py`.
2. Calculate the effective ground temperature using [`calculate_effective_ground_temperature`](@ref).
3. Return the components for the [`WeatherData`](@ref) constructor.
"""
function process_weather(
    archetype::Object,
    scope_data::ScopeData,
    envelope_data::EnvelopeData,
    building_nodes::BuildingNodeNetwork;
    ignore_year::Bool=false,
    repeat::Bool=false,
    save_layouts::Bool=true,
    resampling::Int=5,
    mod::Module=@__MODULE__,
    realization::Symbol=:realization
)
    # Process the weather and preliminary demand data
    heating_demand_W,
    cooling_demand_W,
    ambient_temperature_K,
    total_effective_irradiation_W_effm2,
    heating_set_point_K,
    cooling_set_point_K = create_building_weather(
        archetype,
        scope_data,
        envelope_data,
        building_nodes;
        ignore_year=ignore_year,
        repeat=repeat,
        save_layouts=save_layouts,
        resampling=resampling,
        mod=mod
    )

    # Calculate the effective ground temperature and check it's ok.
    ground_temperature_K = calculate_effective_ground_temperature(ambient_temperature_K)

    # Return the components for WeatherData
    return heating_demand_W,
    cooling_demand_W,
    ambient_temperature_K,
    ground_temperature_K,
    total_effective_irradiation_W_effm2,
    heating_set_point_K,
    cooling_set_point_K
end


"""
    calculate_effective_ground_temperature(
        ambient_temp_K::SpineDataType;
        coeff::Real = 1.7
        realization::Symbol = :realization,
    )

Calculate the effective ground temperature based on the ambient temperature.

The method used in this function is based on 
*Kissock et al., Simplified Model for Ground Heat Transfer from Slab-on-Grade buildings, ASHRAE Transactions 2013*
```math
\\frac{a * T_\\text{ambient annual average} + T_\\text{ambient 3-month moving average}}{1 + a},
```
where the default `coeff` is `a=1.7`.

*NOTE!* The ambient temperature timeseries is assumed to repeat when calculating
the annual and 3-month moving averages.
"""
function calculate_effective_ground_temperature(
    ambient_temp_K::TimeSeries;
    coeff::Real=1.7,
    realization::Symbol=:realization
)
    # Ambient temperature assumed to repeat when calculating moving averages
    repeating_ambient = parameter_value(
        TimeSeries(
            ambient_temp_K.indexes,
            ambient_temp_K.values,
            ambient_temp_K.ignore_year,
            true,
        ),
    )

    # Calculate the moving annual average
    annual_MA_timeslices = [TimeSlice(ts - Year(1), ts) for ts in ambient_temp_K.indexes]
    annual_MA_values = [repeating_ambient(t=ts) for ts in annual_MA_timeslices]
    annual_MA = TimeSeries(
        ambient_temp_K.indexes,
        annual_MA_values,
        ambient_temp_K.ignore_year,
        ambient_temp_K.repeat,
    )

    # Calculate the previous 3-month moving average
    three_month_MA_timeslices =
        [TimeSlice(ts - Month(3), ts) for ts in ambient_temp_K.indexes]
    three_month_MA_values = [repeating_ambient(t=ts) for ts in three_month_MA_timeslices]
    three_month_MA = TimeSeries(
        ambient_temp_K.indexes,
        three_month_MA_values,
        ambient_temp_K.ignore_year,
        ambient_temp_K.repeat,
    )

    # Calculate and return the effective ground temperature
    return (coeff * annual_MA + three_month_MA) / (1 + coeff)
end
function calculate_effective_ground_temperature(
    ambient_temp_K::TimePattern;
    coeff::Real=1.7,
    realization::Symbol=:realization
)
    @error """
    `TimePattern` form ambient temperatures currently unsupported!
    Please use `TimeSeries` instead.
    """
    return nothing
end
function calculate_effective_ground_temperature(
    ambient_temp_K::Map;
    coeff::Real=1.7,
    realization::Symbol=:realization
)
    return calculate_effective_ground_temperature(
        ambient_temp_K[realization];
        coeff=coeff
    )
end
calculate_effective_ground_temperature(
    ambient_temp_K::Real;
    coeff::Real=1.7,
    realization::Symbol=:realization
) = ambient_temp_K


"""
    create_building_weather(
        archetype::Object,
        scope_data::ScopeData,
        envelope_data::EnvelopeData,
        building_nodes::BuildingNodeNetwork;
        ignore_year::Bool=false,
        repeat::Bool=false,
        save_layouts::Bool=true,
        resampling::Int=5,
        mod::Module=@__MODULE__
    )

Try to create `building_weather` automatically using `ArBuWe.py`.

NOTE! The `mod` keyword changes from which Module data is accessed from,
`@__MODULE__` by default.

Essentially, automatically fetches weather data from ERA5 using the
`PYPSA/atlite` python library, and aggregates it according to the GIS data
indicated via the `shapefile_path` and `raster_weight_path` parameters for
the `building_stock` objects. The desired weather period needs to be indicated
using the `building_archetype` `weather_start` and `weather_end` parameters,
and weighting is done based on the `building_scope`, the shapefile at `shapefile_path`,
and the optional raster data at `raster_weight_path`.

The optional `ignore_year` and `repeat` keywords are used to control the
corresponding flags of the created `SpineInterface.TimeSeries`.
By default, the created `TimeSeries` are year-aware and non-repeating.
The `save_layouts` keyword is used to control whether the layouts used for
weighting the weather data are saved for diagnostics.

The actual calculations are handled primarily through the
`aggregate_demand_and_weather` function of `ArBuWe.py`.
"""
function create_building_weather(
    archetype::Object,
    scope_data::ScopeData,
    envelope_data::EnvelopeData,
    building_nodes::BuildingNodeNetwork;
    ignore_year::Bool=false,
    repeat::Bool=false,
    save_layouts::Bool=true,
    resampling::Int=5,
    mod::Module=@__MODULE__
)
    # Fetch air node data, as well as any other nodes with set points.
    (air_node, air_node_data) = only(
        filter(pair -> pair[2].interior_air_and_furniture_weight > 0, building_nodes)
    )
    set_nodes = filter(
        pair -> !isnothing(pair[2].heating_set_point_K), building_nodes
    )
    pop!(set_nodes, air_node)

    # Import `ArBuWe.py`, doesn't work outside the function for some reason...
    abw = pyimport("ArBuWe")
    # Fetch the information necessary for `ArBuWe.py`.
    w_start = string(mod.weather_start(building_archetype=archetype))
    w_end = string(mod.weather_end(building_archetype=archetype))
    weights =
        Dict(string(key.name) => val for (key, val) in scope_data.location_id_gfa_weights)
    bw_name = string(scope_data.building_scope.name) * '_' * w_start * '_' * w_end

    # Convert heating and cooling set points to TimeSeries value arrays for python xarray processing.
    hourly_inds = collect(DateTime(w_start):Hour(1):DateTime(w_end)+Day(31)) # Need to play it safe because of atlite time stamps.
    zero_ts = TimeSeries(hourly_inds, zeros(size(hourly_inds)))
    heating_set_point_K = zero_ts + mod.indoor_air_heating_set_point_override_K(building_archetype=archetype)
    cooling_set_point_K = zero_ts + mod.indoor_air_cooling_set_point_override_K(building_archetype=archetype)

    # Calculate total internal heat gains including connected set point nodes.
    internal_heat_gains_heating_W =
        zero_ts + (
            air_node_data.internal_heat_gains_air_W + sum(
                (
                    set_node_data.heat_transfer_coefficient_structures_interior_W_K +
                    get(set_node_data.heat_transfer_coefficients_base_W_K, air_node, 0.0) +
                    get(set_node_data.heat_transfer_coefficients_gfa_scaled_W_K, air_node, 0.0)
                ) * (
                    set_node_data.heating_set_point_K - heating_set_point_K
                )
                for (set_node, set_node_data) in set_nodes
            )
        )
    internal_heat_gains_cooling_W =
        zero_ts + (
            air_node_data.internal_heat_gains_air_W + sum(
                (
                    set_node_data.heat_transfer_coefficient_structures_interior_W_K +
                    get(set_node_data.heat_transfer_coefficients_base_W_K, air_node, 0.0) +
                    get(set_node_data.heat_transfer_coefficients_gfa_scaled_W_K, air_node, 0.0)
                ) * (
                    set_node_data.cooling_set_point_K - cooling_set_point_K
                )
                for (set_node, set_node_data) in set_nodes
            )
        )

    # Calculate required node properties.
    self_discharge_coefficient_W_K =
        zero_ts +
        air_node_data.self_discharge_base_W_K +
        air_node_data.self_discharge_gfa_scaled_W_K
    total_ambient_heat_transfer_coefficient_with_HRU_W_K =
        zero_ts +
        air_node_data.heat_transfer_coefficient_windows_W_K +
        air_node_data.heat_transfer_coefficient_ventilation_and_infiltration_W_K +
        air_node_data.heat_transfer_coefficient_thermal_bridges_W_K
    total_ambient_heat_transfer_coefficient_without_HRU_W_K =
        zero_ts +
        air_node_data.heat_transfer_coefficient_windows_W_K +
        air_node_data.heat_transfer_coefficient_ventilation_and_infiltration_W_K_HRU_bypass +
        air_node_data.heat_transfer_coefficient_thermal_bridges_W_K

    # Calculate horizontal vs vertical window area.
    horizontal_window_surface_area_m2 =
        envelope_data.window.surface_area_m2 *
        mod.window_area_distribution_towards_cardinal_directions(
            building_archetype=archetype,
            cardinal_direction=:horizontal
        )
    vertical_window_surface_area_m2 =
        envelope_data.window.surface_area_m2 -
        horizontal_window_surface_area_m2

    # Call `ArBuWe.py` to fetch and aggregate the weather.
    heating_demand_W,
    cooling_demand_W,
    ambient_temperature_K,
    total_effective_irradiation_W_effm2 = abw.aggregate_demand_and_weather(
        scope_data.shapefile_path,
        w_start,
        w_end,
        weights,
        mod.external_shading_coefficient(building_archetype=archetype),
        values(heating_set_point_K),
        values(cooling_set_point_K),
        values(internal_heat_gains_heating_W),
        values(internal_heat_gains_cooling_W),
        values(self_discharge_coefficient_W_K),
        values(total_ambient_heat_transfer_coefficient_with_HRU_W_K),
        values(total_ambient_heat_transfer_coefficient_without_HRU_W_K),
        mod.solar_heat_gain_convective_fraction(building_archetype=archetype),
        mod.window_non_perpendicularity_correction_factor(building_archetype=archetype),
        scope_data.total_normal_solar_energy_transmittance,
        vertical_window_surface_area_m2,
        horizontal_window_surface_area_m2,
        scope_data.raster_weight_path,
        resampling,
        bw_name,
        save_layouts,
    )

    # Convert the `PyObjects` to Spine data structures.
    heating_demand_W = _pyseries_to_timeseries(
        heating_demand_W;
        ignore_year=ignore_year,
        repeat=repeat
    )
    cooling_demand_W = _pyseries_to_timeseries(
        cooling_demand_W;
        ignore_year=ignore_year,
        repeat=repeat
    )
    ambient_temperature_K = _pyseries_to_timeseries(
        ambient_temperature_K;
        ignore_year=ignore_year,
        repeat=repeat
    )
    total_effective_irradiation_W_effm2 = Dict(
        Symbol(key) => _pyseries_to_timeseries(
            val;
            ignore_year=ignore_year,
            repeat=repeat
        ) for (key, val) in total_effective_irradiation_W_effm2
    )

    # Match set point time series lengths with weather.
    heating_set_point_K = TimeSeries(
        ambient_temperature_K.indexes,
        heating_set_point_K.values[1:length(ambient_temperature_K.indexes)],
        heating_set_point_K.ignore_year,
        heating_set_point_K.repeat
    )
    cooling_set_point_K = TimeSeries(
        ambient_temperature_K.indexes,
        cooling_set_point_K.values[1:length(ambient_temperature_K.indexes)],
        cooling_set_point_K.ignore_year,
        cooling_set_point_K.repeat
    )
    return heating_demand_W,
    cooling_demand_W,
    ambient_temperature_K,
    total_effective_irradiation_W_effm2,
    heating_set_point_K,
    cooling_set_point_K
end


"""
    _pyseries_to_timeseries(
        pyseries::PyCall.PyObject;
        ignore_year::Bool = false,
        repeat::Bool = false,
    )

Convert `ArBuWe.py` output `pandas.Series` into a `TimeSeries`.

The optional keywords can be used to tweak how the `TimeSeries` is flagged.
Default is a year-aware repeating timeseries.
"""
function _pyseries_to_timeseries(
    pyseries::PyCall.PyObject;
    ignore_year::Bool=false,
    repeat::Bool=false
)
    TimeSeries(collect(pyseries.index), pyseries.values, ignore_year, repeat)
end
