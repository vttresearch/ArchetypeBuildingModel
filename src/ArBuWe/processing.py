# ArBuWe.py

# The main module file, contains the functions for the archetype building weather processing.

from pathlib import Path
import geopandas
import atlite
import atlite.pv.irradiation as irr
import atlite.pv.orientation as ori
import atlite.pv.solar_position as sol
import rioxarray
import numpy as np
import xarray
import matplotlib.pyplot as plt


class Shapefile:
    """A shapefile object for storing the loaded GeoDataFrame."""

    def __init__(self, shapefile_path, resolution=0.01, loose_bound_offset=0.25):
        """
        Create a `Shapefile` object based on an existing shapefile.

        Parameters
        ----------
        shapefile_path : str | path-like
            Filepath pointing to the shapefile.

        Optional parameters
        -------------------
        resolution : float
            Desired resolution of the generated uniform raster, 0.01 degrees by default.
        loose_bound_offset : float
            Raster bounds increased over the shapefile bounds, default 0.25 degrees due to ERA5 resolution.

        """
        self.name = Path(shapefile_path).stem
        self.data = geopandas.read_file(shapefile_path).set_index("location")
        self.resolution = resolution
        self.loose_bound_offset = loose_bound_offset
        self.loose_bounds = self.calculate_loose_bounds()
        self.raster = self.rasterize()

    def calculate_loose_bounds(self):
        """
        Calculates the loose bounds for the shapefile, based on `loose_bound_offset`.

        Returns
        -------
        A tuple with the loosened bounding box of the shapefile.
        """
        # Find and loosen the bounds of the shapefile
        x1, y1, x2, y2 = self.data.total_bounds
        x1 -= self.loose_bound_offset
        x2 += self.loose_bound_offset
        y1 -= self.loose_bound_offset
        y2 += self.loose_bound_offset

        return x1, y1, x2, y2

    def rasterize(self):
        """
        Initialize a uniform raster based on the shapefile.

        Returns
        -------
        raster : An `xarray.DataArray` corresponding to the shapefile.
        """
        # Fetch the loose bounds
        x1, y1, x2, y2 = self.loose_bounds

        # Form the coordinate arrays based on the bounds.
        x = np.arange(x1, x2, self.resolution)
        np.append(x, x[-1] + self.resolution)
        y = np.arange(y1, y2, self.resolution)
        np.append(y, y[-1] + self.resolution)

        # Create the `xarray.DataArray` for the raster and clip it.
        xar = xarray.DataArray(
            np.ones((y.size, x.size)), coords=[y, x], dims=["y", "x"]
        )
        xar.rio.set_crs(self.data.crs)
        xar = xar.rio.clip(self.data.geometry)
        return xar


def prepare_cutout(
    shapefile,
    weather_start,
    weather_end,
    module="era5",
    features=["influx", "temperature"],
):
    """
    Prepares the `atlite` cutout for the weather data.

    Parameters
    ----------
    shapefile : The `Shapefile` used for defining the bounds of the cutout.
    weather_start : String defining the start of the cutout in `yyyy-mm-dd`, month and day can be skipped.
    weather_end: String defining the end of the cutout in `yyyy-mm-dd`, month and day can be skipped.
    module : Module for the weather data, ERA5 by default.
    features : Climate data features to be fetched, `influx` and `temperature` by default.

    Returns
    -------
    cutout : Returns the prepared `atlite` cutout.
    """
    # Load the shapefile and define the bounds for the cutout.
    x1, y1, x2, y2 = shapefile.loose_bounds

    # Define and prepare the cutout.
    cutout = atlite.Cutout(
        path=Path(
            "data/" + shapefile.name + "_" + weather_start + "_" + weather_end
        ).with_suffix(".nc"),
        module=module,
        x=slice(x1, x2),
        y=slice(y1, y2),
        time=slice(weather_start, weather_end),
    )
    cutout.prepare(features=features)

    # Check that cutout and shapefile have the same coordinate reference system
    if shapefile.data.crs != cutout.crs:
        raise ValueError(
            f"""
            The coordinate reference systems of the shapefile and the cutout don't match!
            Shapefile CRS: {shapefile.crs}
            Cutout CRS: {cutout.crs}
            """
        )
    return cutout


def prepare_layout(shapefile, cutout, weights, raster_path=None, resampling=5):
    """
    Prepares a layout for aggregating the weather data.

    Uses `rio.reproject` to resample the underlying raster data to match the
    resolution required by the atlite `cutout`.

    Parameters
    ----------
    shapefile : The shapefile used as the basis for the layout.
    cutout : The atlite cutout to match the layout to.
    weights : A dictionary with weights for the `location_id` fields in the shapefile.

    Other Parameters
    ----------------
    raster_path=None : Optional path to further raster data used for generating the layout.
    resampling=5 : Setting for resampling the data, `5=average` by default.
        See [`rasterio.enums.Resampling`](https://rasterio.readthedocs.io/en/stable/api/rasterio.enums.html#rasterio.enums.Resampling)

    Returns
    -------
    raster : The original raster used for the layout.
    layout : The layout reprojected and resampled to match the atlite `cutout`.
    """
    # Create a reduced GeoDataFrame with only the weighted geometries and include the weight data.
    locations = weights.keys()
    gdf = shapefile.data.loc[locations]
    gdf["weight"] = weights.values()

    # If `raster_path` is defined, load the raster data and clip it with the whole shapefile.
    if raster_path is not None:
        raster = (
            rioxarray.open_rasterio(raster_path, masked=True)
            .rio.clip(shapefile.data.geometry, from_disk=True)
            .squeeze()
        )
    else:  # Else, use the default uniform raster created based on the `Shapefile`
        raster = shapefile.raster

    # Check if raster and cutout CRS matches
    if raster.rio.crs != cutout.crs:
        raise ValueError(
            f"""
            The coordinate reference systems of the cutout and the raster don't match!
            Cutout CRS: {cutout.crs}
            Raster CRS: {raster.rio.crs}
            """
        )

    # Take the given weights into account and normalize the rasters.
    # Reprojection to layout resolution is done on sub-raster basis to avoid memory issues.
    rasters = []
    for lid in locations:
        # Define a loose bounding box based on shapefile settings
        minx, miny, maxx, maxy = gdf.geometry.loc[lid].bounds
        minx -= shapefile.loose_bound_offset
        miny -= shapefile.loose_bound_offset
        maxx += shapefile.loose_bound_offset
        maxy += shapefile.loose_bound_offset

        # Clip the loose box from the original raster and drop the rest for computational efficiency
        rst = raster.rio.clip_box(minx, miny, maxx, maxy)
        # Clip the municipality tightly, but don't `drop` to avoid issues with resampling.
        rst = (
            rst.rio.clip(
                [gdf.geometry.loc[lid]], all_touched=True, drop=False, from_disk=True
            )
            .fillna(0.0)
            .rio.reproject(
                cutout.crs,
                shape=cutout.shape,
                transform=cutout.transform,
                resampling=resampling,
                from_disk=True,
                nodata=0.0,
            )
        )
        rst = rst / rst.sum() * weights[lid]  # Normalize and weigh the sub-shape.
        rasters.append(rst)  # Collect the weighted sub-raster into the list.

    # Combine the weighted and normalized rasters into one.
    layout = sum(rasters)
    layout = layout / layout.sum()  # Finally, normalize the whole layout

    return raster, layout


##TODO: REMOVE OLD CODE!
# def process_weather(cutout, layout):
#     """
#     Processes the weather data using the given layout.

#     The intended use for this function is to calculate the weighted average
#     weather parameters for the `ArBuMo`, but depending on
#     the used layout it can be used to do other things as well.

#     Parameters
#     ----------
#     cutout : `atlite` cutout containing the weather data.
#     layout : `xarray.DataArray` for the desired layout.

#     Returns
#     -------
#     ambient_temperature : Ambient temperature timeseries [K].
#     diffuse_irradiation : Diffuse irradiation timeseries [W/m2].
#     direct_irradiation :
#         Direct irradiation on walls facing the cardinal directions [W/m2].
#     """
#     # Define the slope and azimuth angles for the horizontal and cardinal directions.
#     dirs = {
#         "horizontal": (0.0, 0.0),
#         "north": (90.0, 0.0),
#         "east": (90.0, 90.0),
#         "south": (90.0, 180.0),
#         "west": (90.0, 270.0),
#     }

#     # Calculate the aggregated weather using the `atlite` cutout.
#     ambient_temperature = (
#         cutout.temperature(layout=layout).squeeze().to_series() + 273.15
#     )
#     diffuse_irradiation = (
#         cutout.irradiation(
#             orientation={"slope": 0.0, "azimuth": 0.0},
#             layout=layout,
#             irradiation="diffuse",
#         )
#         .squeeze()
#         .to_series()
#     )
#     direct_irradiation = {
#         dir: cutout.irradiation(
#             orientation={"slope": sl, "azimuth": az},
#             layout=layout,
#             irradiation="direct",
#         )
#         .squeeze()
#         .to_series()
#         for dir, (sl, az) in dirs.items()
#     }

#     return ambient_temperature, diffuse_irradiation, direct_irradiation


def plot_layout(shapefile, layout, dpi=300, title="layout"):
    """
    Plot layout to see if it makes sense.

    Parameters
    ----------
    shapefile : Shapefile
        Shapefile used for creating the layout.
    layout : xarray.DataArray
        The layout to be plotted

    Optional parameters
    -------------------
    dpi : int
        Dots-per-inch for the output figure
    title : str
        Title for the figure

    Returns
    -------
    A pyplot figure.
    """
    f, ax = plt.subplots(dpi=dpi)
    try:
        layout.where(layout != xarray.zeros_like(layout)).plot(ax=ax)
    except:
        print("Cannot print point layout raster!")
    shapefile.data.geometry.plot(
        ax=ax, edgecolor=(0, 0, 0, 1), facecolor=(0, 0, 0, 0), linewidth=0.2
    )
    ax.set_title(title)
    return f


## TODO: REMOVE OLD CODE!
# def aggregate_weather(
#     shapefile_path,
#     weights,
#     weather_start,
#     weather_end,
#     raster_path=None,
#     filename="scope",
#     save_layouts=True,
#     resampling=5,
# ):
#     """
#     Aggregate weather data for the given input arguments.

#     The shapefile is used to define the geographical area under investigation,
#     as well as to assign the `weights` to its different areas.
#     If the `raster_path` is included, its data will be included when generating
#     the layout for weighting the weather data.
#     Otherwise, an uniform base distribution is assumed for the different
#     areas in the shapefile.

#     Parameters
#     ----------
#     shapefile_path : str | path-like
#         Path to the shapefile defiining the area under investigation.
#     weights : dict
#         Weights for the different polygons in the shapefile.
#     weather_start : str | datetime-like
#         Time when weather data period starts in `yyyy-mm-dd`, month and day can be omitted.
#     weather_end : str | datetime-like
#         Time when weather data period ends in `yyyy-mm-dd`, month and day can be omitted.

#     Optional parameters
#     -------------------
#     raster_path = None : str | path-like
#         Path to an optional raster file to include in the weighting.
#     filename = "scope" : str
#         Name of the image file displaying the final layout.
#     save_layouts = True : bool
#         A flag to control whether layout images are to be saved.
#     resampling=5 : Setting for resampling the data for `prepare_layout()`, `5=average` by default.

#     Returns
#     -------
#     ambient_temperature : Ambient temperature timeseries [K].
#     diffuse_irradiation : Diffuse irradiation timeseries [W/m2].
#     direct_irradiation :
#         Direct irradiation on walls facing the cardinal directions [W/m2].
#     """
#     shapefile = Shapefile(shapefile_path)
#     cutout = prepare_cutout(shapefile, weather_start, weather_end)
#     raster, layout = prepare_layout(shapefile, cutout, weights, raster_path, resampling)
#     if save_layouts:
#         for rst, name in [(raster, "raster"), (layout, "layout")]:
#             f = plot_layout(shapefile, rst, dpi=600, title=filename + " " + name)
#             f.savefig("figs/" + filename + "_" + name + ".png", dpi=600)
#             plt.close(f)
#     return process_weather(cutout, layout)


def preprocess_weather(cutout, external_shading_coefficient):
    """
    Preprocesses weather data for heating and cooling demand calculations.

    This function calculates the aggregated weather quantities
    required for the heating/cooling demand calculations.
    Note that the returned `effective_irradiation` contains the
    effective total irradiation for horizontal and vertical surfaces
    respectively. These are direction-agnostic, and designed to be used
    with the full horizontal/vertical envelope surface areas respectively.
    For the vertical surfaces, we're assuming that
    they are distributed equally to all directions *(like a circle)*.
    Thus, the effective direct irradiation on the equidistributed
    surface is calculated as `I_dir_eff = c_shading * I_dir / Pi`.

    Parameters
    ----------
    cutout : `atlite` cutout containing the weather data.
    external_shading_coefficient : `float`
        Coefficient for direct irradiation quantities to account for shading due to environment.

    Returns
    -------
    ambient_temperature : `xarray.DataSet`
        Ambient temperature timeseries [K].
    effective_irradiation : `dict`
        Total effective irradiation timeseries [W/m2] for horizontal and vertical surfaces.
    """
    # Define the slope and azimuth angles for the horizontal and vertical surfaces.
    dirs = {  # Azimuths don't really matter due to slope and tracking respectively.
        "horizontal": ori.get_orientation({"slope": 0.0, "azimuth": 0.0}),
        "vertical": ori.get_orientation({"slope": 90.0, "azimuth": 0.0}),
    }

    # Get solar position and surface orientations
    solar_position = sol.SolarPosition(cutout.data)
    surface_orientation = {
        dir: ori.SurfaceOrientation(
            cutout.data,
            solar_position,
            orientation,
            "vertical",
        )
        for dir, orientation in dirs.items()
    }

    # Calculate total diffuse irradiation (diffuse + ground) and direct irradiation.
    diffuse_irradiation_W_m2 = {
        dir: irr.TiltedIrradiation(
            cutout.data,
            solar_position,
            orientation,
            trigon_model="simple",
            clearsky_model="simple",
            tracking="vertical",
            irradiation="diffuse",
        )
        + irr.TiltedIrradiation(
            cutout.data,
            solar_position,
            orientation,
            trigon_model="simple",
            clearsky_model="simple",
            tracking="vertical",
            irradiation="ground",
        )
        for dir, orientation in surface_orientation.items()
    }
    direct_irradiation_W_m2 = {
        dir: irr.TiltedIrradiation(
            cutout.data,
            solar_position,
            orientation,
            trigon_model="simple",
            clearsky_model="simple",
            tracking="vertical",
            irradiation="direct",
        )
        for dir, orientation in surface_orientation.items()
    }

    # Calculate the total effective irradiations for horizontal and vertical surfaces respectively.
    # Note that the `external_shading_coefficient` only applies to direct irradiation!
    # Also note that direct irradiation is divided per `np.pi` to account for not all vertical surface area facing the right way!
    total_effective_irradiation_W_effm2 = {
        dir: diff + external_shading_coefficient * direct_irradiation_W_m2[dir] / np.pi
        for dir, diff in diffuse_irradiation_W_m2.items()
    }

    # Ambient temperatures
    ambient_temperature_K = cutout.data["temperature"]

    # Return the weather quantities
    return ambient_temperature_K, total_effective_irradiation_W_effm2


def expand_to_xarray(array, xarray_to_match, description, units):
    """
    Expands an array of values to match a given xarray.

    Essentially clones the time series data in the given array
    into an xarray with dimensions matching the given xarray.

    Parameters
    ----------
    array : np.array
        The data to be expanded into an xarray.
    xarray_to_match : xarray.DataArray
        The data to match when expanding.
    description : string
        A brief description of the data to be expanded.
    units : string
        The unit of the data values.

    Returns
    -------
    xarray : xarray.DataArray
        The original array extended to a matching xarray.
    """
    # Fetch dimensions from the xarray to match
    lon = xarray_to_match.x
    lat = xarray_to_match.y
    time = xarray_to_match.time

    # Check that the time dimension lengths match.
    if len(array) != len(time):
        raise ValueError(
            f"""
            The temporal lengths of `array` and `xarray_to_match` don't match!
            Array length: {len(array)}
            Xarray time length: {len(time)}
            """
        )

    # Expand data to the desired size and create the xarray
    data = np.repeat(array, len(lon) * len(lat)).reshape(xarray_to_match.shape)
    xar = xarray.DataArray(
        data=data,
        dims=xarray_to_match.dims,
        coords=xarray_to_match.coords,
        attrs=dict(description=description, units=units),
    )
    return xar


def process_initial_heating_demand(
    set_point_K,
    ambient_temperature_K,
    total_effective_irradiation_W_effm2,
    internal_heat_gains_W,
    self_discharge_coefficient_W_K,
    total_ambient_heat_transfer_coefficient_W_K,
    solar_heat_gain_convective_fraction,
    window_non_perpendicularity_correction_factor,
    total_normal_solar_energy_transmittance,
    vertical_window_surface_area_m2,
    horizontal_window_surface_area_m2,
):
    """
    Processes the preliminary aggregated heating demand.

    TODO: DOCUMENTATION!

    Parameters
    ----------

    Returns
    -------

    """
    initial_heating_demand = (
        self_discharge_coefficient_W_K * set_point_K  # Self-discharge heat losses
        + total_ambient_heat_transfer_coefficient_W_K  # Ambient heat losses: This need to account for HRU bypass for cooling!
        * (set_point_K - ambient_temperature_K)  # Ambient heat losses.
        - internal_heat_gains_W  # NOTE! These need to include utilisable losses from the DHW tank!
        - solar_heat_gain_convective_fraction  # Solar gains are more complicated.
        * window_non_perpendicularity_correction_factor
        * total_normal_solar_energy_transmittance
        * (
            vertical_window_surface_area_m2
            * total_effective_irradiation_W_effm2["vertical"]
            + horizontal_window_surface_area_m2
            * total_effective_irradiation_W_effm2["horizontal"]
        )
    )
    return initial_heating_demand
