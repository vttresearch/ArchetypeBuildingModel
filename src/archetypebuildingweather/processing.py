# archetypebuildingweather.py

# The main module file, contains the functions for the archetype building weather processing.

from pathlib import Path
import geopandas
import atlite
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
    shapefile, weather_start, weather_end, module="era5", features=["influx", "temperature"]
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
        path=Path("data/" + shapefile.name + "_" + weather_start + "_" + weather_end).with_suffix(".nc"),
        module=module,
        x=slice(x1, x2),
        y=slice(y1, y2),
        time=slice(weather_start, weather_end),
    )
    cutout.prepare(features=features)
    return cutout


def coarsen_raster(raster, target, method="mean", boundary="trim"):
    """
    Coarsens `raster` data to match the resolution of the `target` raster.

    Parameters
    ----------
    raster : The `xarray.DataArray` to be coarsened.
    target : The `xarray.DataArray` to match.

    Optional parameters
    -------------------
    method : Default 'mean' for averaging, 'sum' for summation.
    boundary : Setting for `xarray.coarsen()` sampling, `trim` by default.
    """
    # First, figure out how much we need to coarsen the data.
    if raster.x.values.size > 1:
        x_c = int(abs(np.diff(target.x.values).mean() // np.diff(raster.x.values).mean()))
    else:
        x_c = 1
    if raster.y.values.size > 1:
        y_c = int(abs(np.diff(target.y.values).mean() // np.diff(raster.y.values).mean()))
    else:
        y_c = 1

    # Next, coarsen the raster.
    if method == "mean":
        coarse_raster = raster.coarsen(x=x_c, y=y_c, boundary=boundary).mean()
    elif method == "sum":
        coarse_raster = raster.coarsen(x=x_c, y=y_c, boundary=boundary).sum()
    else:
        raise "`method` not recognized! Please use either 'mean' or 'sum'."

    return coarse_raster


def prepare_layout(
    shapefile, weights, raster_path=None, method="mean", boundary="trim"
):
    """
    Prepares a layout for aggregating the weather data.

    Note that the layout generated by this function uses the resolution of the
    underlying raster data, and isn't compatible with `atlite` as is.
    See `match_layout` for resampling the data to match the `atlite` cutout.

    Parameters
    ----------
    shapefile : The shapefile used as the basis for the layout.
    weights : A dictionary with weights for the `location_id` fields in the shapefile.

    Other Parameters
    ----------------
    raster_path : Optional path to further raster data used for generating the layout.
    method : Default 'mean' for averaging, 'sum' for summation of raster coarsening.
    boundary : Setting for `xarray.coarsen()` sampling, `trim` by default.
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
        # Coarsen the raster to conserve memory.
        raster = coarsen_raster(
            raster, shapefile.raster, method=method, boundary=boundary
        )

    else:  # Else, use the default uniform raster created based on the `Shapefile`
        raster = shapefile.raster

    # Take the given weights into account and normalize the rasters.
    rasters = []
    for lid in locations:
        rst = raster.rio.clip(
            [gdf.geometry.loc[lid]], all_touched=True, drop=False, from_disk=True
        ).fillna(0.0)
        rst = rst / rst.sum() * weights[lid]
        rasters.append(rst)

    # Combine the weighted and normalized rasters into one.
    raster = sum(rasters)
    return raster


def match_layout(
    shapefile,
    raster,
    cutout,
    boundary="trim",
    coarse_method="mean",
    reindex_method="nearest",
):
    """
    Resamples and normalizes a layout to match the `atlite` cutout.

    Parameters
    ----------
    shapefile : The original shapefile, used for clipping to mitigate grid mismatch.
    raster : The desired `xarray.DataArray` layout to match.
    cutout : The `atlite.Cutout` to match the `raster` to.

    Optional parameters
    -------------------
    boundary : Setting for `xarray.coarsen()` sampling, `trim` by default.
    coarse_method : Default 'mean' for averaging, 'sum' for summation of raster coarsening.
    reindex_method : Setting for `xarray.reindex_like()`, `nearest` by default.

    Returns
    -------
    layout : The resampled and normalized layout.
    """
    # Coarsen and reindex the raster.
    layout = coarsen_raster(
        raster, cutout.data, method=coarse_method, boundary=boundary
    ).reindex_like(cutout.data, method=reindex_method)

    # Clip the reindexed raster to mitigate grid mismatch.
    layout = layout.rio.clip(
        shapefile.data.geometry, all_touched=True, drop=False, from_disk=True
    ).fillna(0.0)

    # Finally, normalize the layout
    layout = layout / layout.sum()

    return layout


def process_weather(cutout, layout):
    """
    Processes the weather data using the given layout.

    The intended use for this function is to calculate the weighted average
    weather parameters for the `ArchetypeBuildingModel`, but depending on
    the used layout it can be used to do other things as well.

    Parameters
    ----------
    cutout : `atlite` cutout containing the weather data.
    layout : `xarray.DataArray` for the desired layout.

    Returns
    -------
    ambient_temperature : Ambient temperature timeseries [K].
    diffuse_irradiation : Diffuse irradiation timeseries [W/m2].
    direct_irradiation :
        Direct irradiation on walls facing the cardinal directions [W/m2].
    """
    # Define the slope and azimuth angles for the horizontal and cardinal directions.
    dirs = {
        "horizontal": (0.0, 0.0),
        "north": (90.0, 0.0),
        "east": (90.0, 90.0),
        "south": (90.0, 180.0),
        "west": (90.0, 270.0),
    }

    # Calculate the aggregated weather using the `atlite` cutout.
    ambient_temperature = (
        cutout.temperature(layout=layout).squeeze().to_series() + 273.15
    )
    diffuse_irradiation = (
        cutout.irradiation(
            orientation={"slope": 0.0, "azimuth": 0.0},
            layout=layout,
            irradiation="diffuse",
        )
        .squeeze()
        .to_series()
    )
    direct_irradiation = {
        dir: cutout.irradiation(
            orientation={"slope": sl, "azimuth": az},
            layout=layout,
            irradiation="direct",
        )
        .squeeze()
        .to_series()
        for dir, (sl, az) in dirs.items()
    }

    return ambient_temperature, diffuse_irradiation, direct_irradiation


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


def aggregate_weather(
    shapefile_path,
    weights,
    weather_start,
    weather_end,
    raster_path=None,
    filename="scope",
    save_layouts=True,
):
    """
    Aggregate weather data for the given input arguments.

    The shapefile is used to define the geographical area under investigation,
    as well as to assign the `weights` to its different areas.
    If the `raster_path` is included, its data will be included when generating
    the layout for weighting the weather data.
    Otherwise, an uniform base distribution is assumed for the different
    areas in the shapefile.

    Parameters
    ----------
    shapefile_path : str | path-like
        Path to the shapefile defiining the area under investigation.
    weights : dict
        Weights for the different polygons in the shapefile.
    weather_start : str | datetime-like
        Time when weather data period starts in `yyyy-mm-dd`, month and day can be omitted.
    weather_end : str | datetime-like
        Time when weather data period ends in `yyyy-mm-dd`, month and day can be omitted.

    Optional parameters
    -------------------
    raster_path = None : str | path-like
        Path to an optional raster file to include in the weighting.
    filename = "scope" : str
        Name of the image file displaying the final layout.
    save_layouts = True : bool
        A flag to control whether layout images are to be saved.

    Returns
    -------
    ambient_temperature : Ambient temperature timeseries [K].
    diffuse_irradiation : Diffuse irradiation timeseries [W/m2].
    direct_irradiation :
        Direct irradiation on walls facing the cardinal directions [W/m2].
    """
    shapefile = Shapefile(shapefile_path)
    cutout = prepare_cutout(shapefile, weather_start, weather_end)
    raster = prepare_layout(shapefile, weights, raster_path)
    layout = match_layout(shapefile, raster, cutout)
    if save_layouts:
        plot_layout(shapefile, raster, dpi=600, title=filename + " raster").savefig(
            "figs/" + filename + "_raster.png", dpi=600
        )
        plot_layout(shapefile, layout, dpi=600, title=filename + " layout").savefig(
            "figs/" + filename + "_layout.png", dpi=600
        )
    return process_weather(cutout, layout)
