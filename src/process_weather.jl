#=
    process_weather.jl

Contains functions for processing weather data.
=#


"""
    process_weather(weather::Object; mod::Module = @__MODULE__)

Process `weather` data for the [`WeatherData`](@ref) constructor.

NOTE! The `mod` keyword changes from which Module data is accessed from,
`@__MODULE__` by default.

Essentially, performs the following steps:
1. Fetch the ambient temperature data for `weather`.
2. Calculate the effective ground temperature using [`calculate_effective_ground_temperature`](@ref).
3. Fetch diffuse and direct solar irradiation data.
4. Return the components for the [`WeatherData`](@ref) constructor.
"""
function process_weather(weather::Object; mod::Module = @__MODULE__)
    # Fetch ambient temperature data and check that it's ok.
    ambient_temp_K = mod.ambient_temperature_K(building_weather = weather)
    all(values(ambient_temp_K) .>= 0) || @warn """
    `ambient_temperature_K` for `$(weather)` shouldn't have negative values!
    $(count(values(ambient_temp_K) .< 0)) violations found,
    with a minimum of $(minimum(values(ambient_temp_K))).
    """

    # Calculate the effective ground temperature and check it's ok.
    ground_temp_K = calculate_effective_ground_temperature(ambient_temp_K)
    all(values(ground_temp_K) .>= 0) || @warn """
    `ground_temperature_K` for `$(weather)` shouldn't have negative values!
    $(count(values(ground_temp_K) .< 0)) violations found,
    with a minimum of $(minimum(values(ground_temp_K))).
    """

    # Fetch diffuse solar irradiation data and check it's ok.
    diff_solar_irradiation_W_m2 =
        mod.diffuse_solar_irradiation_W_m2(building_weather = weather)
    all(values(diff_solar_irradiation_W_m2) .>= 0) || @warn """
    `diffuse_solar_irradiation_W_m2` for `$(weather)` shouldn't have negative values!
    $(count(values(diff_solar_irradiation_W_m2) .< 0)) violations found,
    with a minimum of $(minimum(values(diff_solar_irradiation_W_m2))).
    """

    # Fetch direct solar irradiation data and check it's ok.
    dir_solar_irradiation_W_m2 = Dict(
        dir => mod.direct_solar_irradiation_W_m2(
            building_weather = weather;
            cardinal_direction = dir,
        ) for dir in solar_directions
    )
    for dir in solar_directions
        all(values(dir_solar_irradiation_W_m2[dir]) .>= 0) || @warn """
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
        coeff::Real=1.7
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
    coeff::Real = 1.7,
)
    # Ambient temperature assumed to repeat when calculating moving averages
    repeating_ambient = TimeSeries(
        ambient_temp_K.indexes,
        ambient_temp_K.values,
        ambient_temp_K.ignore_year,
        true,
    )

    # Calculate the moving annual average
    annual_MA_timeslices = [TimeSlice(ts - Year(1), ts) for ts in ambient_temp_K.indexes]
    annual_MA_values =
        [parameter_value(repeating_ambient)(ts) for ts in annual_MA_timeslices]
    annual_MA = TimeSeries(
        ambient_temp_K.indexes,
        annual_MA_values,
        ambient_temp_K.ignore_year,
        ambient_temp_K.repeat,
    )

    # Calculate the previous 3-month moving average
    three_month_MA_timeslices =
        [TimeSlice(ts - Month(3), ts) for ts in ambient_temp_K.indexes]
    three_month_MA_values =
        [parameter_value(repeating_ambient)(ts) for ts in three_month_MA_timeslices]
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
    coeff::Real = 1.7,
)
    @error """
    `TimePattern` form ambient temperatures currently unsupported!
    Please use `TimeSeries` instead.
    """
    return nothing
end
function calculate_effective_ground_temperature(ambient_temp_K::Map; coeff::Real = 1.7)
    @error """
    `Map` form ambient temperatures currently unsupported!
    Please use `TimeSeries` instead.
    """
    return nothing
end
calculate_effective_ground_temperature(ambient_temp_K::Real; coeff::Real = 1.7) =
    ambient_temp_K


"""
    create_building_weather(
        archetype::Object,
        scopedata::ScopeData;
        ignore_year::Bool = false,
        repeat::Bool = true,
        save_layouts::Bool = true,
        mod::Module = @__MODULE__,
    )

Try to create `building_weather` automatically using `ArchetypeBuildingWeather.py`.

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
    scopedata::ScopeData;
    ignore_year::Bool = false,
    repeat::Bool = true,
    save_layouts::Bool = true,
    mod::Module = @__MODULE__,
)
    # Import `ArchetypeBuildingWeather.py`, doesn't work outside the function for some reason...
    abw = pyimport("archetypebuildingweather")
    # Fetch the information necessary for `ArchetypeBuildingWeather.py`.
    w_start = string(mod.weather_start(building_archetype = archetype))
    w_end = string(mod.weather_end(building_archetype = archetype))
    weights =
        Dict(string(key.name) => val for (key, val) in scopedata.location_id_gfa_weights)
    bw_name = string(scopedata.building_scope.name) * '_' * w_start * '_' * w_end

    # Call `ArchetypeBuildingWeather.py` to fetch and aggregate the weather.
    ambient_temperature, diffuse_irradiation, direct_irradiation = abw.aggregate_weather(
        scopedata.shapefile_path,
        weights,
        w_start,
        w_end,
        scopedata.raster_weight_path,
        bw_name,
        save_layouts,
    )

    # Convert the `PyObjects` to Spine data structures.
    ambient_temperature = _pyseries_to_timeseries(
        ambient_temperature;
        ignore_year = ignore_year,
        repeat = repeat,
    )
    diffuse_irradiation = _pyseries_to_timeseries(
        diffuse_irradiation;
        ignore_year = ignore_year,
        repeat = repeat,
    )
    direct_irradiation = Map(
        Symbol.(keys(direct_irradiation)),
        _pyseries_to_timeseries.(
            values(direct_irradiation);
            ignore_year = ignore_year,
            repeat = repeat,
        ),
    )

    # Create a new `building_weather` object and its parameter value dictionary
    bw = Object(Symbol(bw_name), :building_weather)
    bw_param_dict = Dict(
        bw => Dict(
            :ambient_temperature_K => parameter_value(ambient_temperature),
            :diffuse_solar_irradiation_W_m2 => parameter_value(diffuse_irradiation),
            :direct_solar_irradiation_W_m2 => parameter_value(direct_irradiation),
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

Convert `ArchetypeBuildingWeather.py` output `pandas.Series` into a `TimeSeries`.

The optional keywords can be used to tweak how the `TimeSeries` is flagged.
Default is a year-aware repeating timeseries.
"""
function _pyseries_to_timeseries(
    pyseries::PyCall.PyObject;
    ignore_year::Bool = false,
    repeat::Bool = true,
)
    TimeSeries(collect(pyseries.index), pyseries.values, ignore_year, repeat)
end