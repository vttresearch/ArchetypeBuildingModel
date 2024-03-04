#=
    process_weather.jl

Contains functions for processing weather data.
=#


"""
    process_weather(
        weather::Object;
        mod::Module = @__MODULE__,
        realization::Symbol = realization,
    )

Process `weather` data for the [`WeatherData`](@ref) constructor.

TODO: Revise documentation!

NOTE! The `mod` keyword changes from which Module data is accessed from,
`@__MODULE__` by default. The `realization` scenario is required for processing
of effective ground temperature.

Essentially, performs the following steps:
1. Fetch the ambient temperature data for `weather`.
2. Calculate the effective ground temperature using [`calculate_effective_ground_temperature`](@ref).
3. Fetch diffuse and direct solar irradiation data.
4. Return the components for the [`WeatherData`](@ref) constructor.
"""
function process_weather(
    weather::Object;
    mod::Module=@__MODULE__,
    realization::Symbol=realization
)
    # Fetch ambient temperature data and check that it's ok.
    ambient_temp_K = mod.ambient_temperature_K(building_weather=weather)
    all(collect_leaf_values(ambient_temp_K) .>= 0) || @warn """
    `ambient_temperature_K` for `$(weather)` shouldn't have negative values!
    $(count(collect_leaf_values(ambient_temp_K) .< 0)) violations found,
    with a minimum of $(minimum(collect_leaf_values(ambient_temp_K))).
    """

    # Calculate the effective ground temperature and check it's ok.
    ground_temp_K = calculate_effective_ground_temperature(ambient_temp_K)
    all(collect_leaf_values(ground_temp_K) .>= 0) || @warn """
    `ground_temperature_K` for `$(weather)` shouldn't have negative values!
    $(count(collect_leaf_values(ground_temp_K) .< 0)) violations found,
    with a minimum of $(minimum(collect_leaf_values(ground_temp_K))).
    """

    # Fetch diffuse solar irradiation data and check it's ok.
    diff_solar_irradiation_W_m2 =
        mod.diffuse_solar_irradiation_W_m2(building_weather=weather)
    all(collect_leaf_values(diff_solar_irradiation_W_m2) .>= 0) || @warn """
    `diffuse_solar_irradiation_W_m2` for `$(weather)` shouldn't have negative values!
    $(count(collect_leaf_values(diff_solar_irradiation_W_m2) .< 0)) violations found,
    with a minimum of $(minimum(collect_leaf_values(diff_solar_irradiation_W_m2))).
    """

    # Fetch direct solar irradiation data and check it's ok.
    dir_solar_irradiation_W_m2 = Dict(
        dir => mod.direct_solar_irradiation_W_m2(
            building_weather=weather;
            cardinal_direction=dir
        ) for dir in solar_directions
    )
    for dir in solar_directions
        all(collect_leaf_values(dir_solar_irradiation_W_m2[dir]) .>= 0) || @warn """
        `direct_solar_irradiation_W_m2[$(dir)]` for `$(weather)` shouldn't have negative values!
        $(count(values(dir_solar_irradiation_W_m2[dir]) .< 0)) violations found,
        with a minimum of $(minimum(values(dir_solar_irradiation_W_m2[dir]))).
        """
    end

    # Return the components for WeatherData
    return ambient_temp_K,
    ground_temp_K,
    diff_solar_irradiation_W_m2,
    dir_solar_irradiation_W_m2
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
        scope_data::ScopeData;
        ignore_year::Bool = false,
        repeat::Bool = true,
        save_layouts::Bool = true,
        resampling::Int=5,
        mod::Module = @__MODULE__,
    )

Try to create `building_weather` automatically using `ArBuWe.py`.

TODO: Revise documentation!

NOTE! The `mod` keyword changes from which Module data is accessed from,
`@__MODULE__` by default.

Essentially tries to automatically fetch weather data from ERA5 using the
`PYPSA/atlite` python library, and aggregate it according to the GIS data
indicated via the `shapefile_path` and `raster_weight_path` parameters for
the `building_stock` objects. The desired weather period needs to be indicated
using the `building_archetype` `weather_start` and `weather_end` parameters,
and weighting is done based on the `building_scope`, the shapefile at `shapefile_path`,
and the optional raster data at `raster_weight_path`.

The optional `ignore_year` and `repeat` keywords are used to control the
corresponding flags of the created `SpineInterface.TimeSeries`.
By default, the created `TimeSeries` are year-aware and repeating.
The `save_layouts` keyword is used to control whether the layouts used for
weighting the weather data are saved for diagnostics.

Returns a new `building_weather` object, as well as a dictionary containing
its parameter values.
"""
function create_building_weather(
    archetype::Object,
    scope_data::ScopeData,
    envelope_data::EnvelopeData,
    building_nodes::BuildingNodeNetwork,
    loads_data::LoadsData;
    ignore_year::Bool=false,
    repeat::Bool=true,
    save_layouts::Bool=true,
    resampling::Int=5,
    mod::Module=@__MODULE__
)
    # Import `ArBuWe.py`, doesn't work outside the function for some reason...
    abw = pyimport("ArBuWe")
    # Fetch the information necessary for `ArBuWe.py`.
    w_start = string(mod.weather_start(building_archetype=archetype))
    w_end = string(mod.weather_end(building_archetype=archetype))
    weights =
        Dict(string(key.name) => val for (key, val) in scope_data.location_id_gfa_weights)
    bw_name = string(scope_data.building_scope.name) * '_' * w_start * '_' * w_end

    # Identify the indoor air and dhw nodes.
    # TODO: Remove the weights, replace with booleans?
    air_node = only(filter(n -> n.interior_air_and_furniture_weight > 0, collect(values(building_nodes))))
    dhw_node = only(filter(n -> n.domestic_hot_water_demand_W != 0, collect(values(building_nodes))))

    # Convert heating and cooling set points to TimeSeries value arrays for python xarray processing.
    hourly_inds = collect(DateTime(w_start):Hour(1):DateTime(w_end)+Day(31)) # Need to play it safe because of atlite time stamps.
    heating_set_point_K = TimeSeries(hourly_inds, zeros(size(hourly_inds))) + mod.indoor_air_heating_set_point_override_K(building_archetype=archetype)
    cooling_set_point_K = TimeSeries(hourly_inds, zeros(size(hourly_inds))) + mod.indoor_air_cooling_set_point_override_K(building_archetype=archetype)

    # Calculate total internal heat gains including utilisable DHW losses.
    # DHW tank heat losses are slightly different for heating and cooling set points...
    internal_heat_gains_heating_W = (
        loads_data.internal_heat_gains_W +
        (dhw_node.minimum_temperature_K - heating_set_point_K) * (
            dhw_node.heat_transfer_coefficients_base_W_K[air_node.building_node] +
            dhw_node.heat_transfer_coefficients_gfa_scaled_W_K[air_node.building_node]
        )
    )
    internal_heat_gains_cooling_W = (
        loads_data.internal_heat_gains_W +
        (dhw_node.minimum_temperature_K - cooling_set_point_K) * (
            dhw_node.heat_transfer_coefficients_base_W_K[air_node.building_node] +
            dhw_node.heat_transfer_coefficients_gfa_scaled_W_K[air_node.building_node]
        )
    )

    # Calculate required node properties.
    self_discharge_coefficient_W_K =
        air_node.self_discharge_base_W_K + air_node.self_discharge_gfa_scaled_W_K
    total_ambient_heat_transfer_coefficient_with_HRU_W_K =
        air_node.heat_transfer_coefficient_windows_W_K +
        air_node.heat_transfer_coefficient_ventilation_and_infiltration_W_K +
        air_node.heat_transfer_coefficient_thermal_bridges_W_K
    total_ambient_heat_transfer_coefficient_without_HRU_W_K =
        air_node.heat_transfer_coefficient_windows_W_K +
        air_node.heat_transfer_coefficient_ventilation_and_infiltration_W_K_HRU_bypass +
        air_node.heat_transfer_coefficient_thermal_bridges_W_K

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
        self_discharge_coefficient_W_K,
        total_ambient_heat_transfer_coefficient_with_HRU_W_K,
        total_ambient_heat_transfer_coefficient_without_HRU_W_K,
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
    total_effective_irradiation_W_effm2 = Map(
        Symbol.(keys(total_effective_irradiation_W_effm2)),
        _pyseries_to_timeseries.(
            values(total_effective_irradiation_W_effm2);
            ignore_year=ignore_year,
            repeat=repeat
        ),
    )

    # Create a new `building_weather` object and its parameter value dictionary
    bw = Object(Symbol(bw_name), :building_weather)
    bw_param_dict = Dict(
        bw => Dict(
            :heating_demand_W => parameter_value(heating_demand_W),
            :cooling_demand_W => parameter_value(cooling_demand_W),
            :ambient_temperature_K => parameter_value(ambient_temperature_K),
            :total_effective_irradiation_W_effm2 => parameter_value(total_effective_irradiation_W_effm2)
        ),
    )
    return bw, bw_param_dict
end


"""
    _pyseries_to_timeseries(
        pyseries::PyCall.PyObject;
        ignore_year::Bool = false,
        repeat::Bool = true,
    )

Convert `ArBuWe.py` output `pandas.Series` into a `TimeSeries`.

The optional keywords can be used to tweak how the `TimeSeries` is flagged.
Default is a year-aware repeating timeseries.
"""
function _pyseries_to_timeseries(
    pyseries::PyCall.PyObject;
    ignore_year::Bool=false,
    repeat::Bool=true
)
    TimeSeries(collect(pyseries.index), pyseries.values, ignore_year, repeat)
end
