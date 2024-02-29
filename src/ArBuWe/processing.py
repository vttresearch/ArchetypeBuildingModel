# ArBuWe.py

# The main module file, contains the functions for the archetype building weather processing.

from pathlib import Path
from geopandas import read_file
from atlite import Cutout
from atlite.pv.irradiation import TiltedIrradiation
from atlite.pv.orientation import SurfaceOrientation, get_orientation
from atlite.pv.solar_position import SolarPosition
from atlite.aggregate import aggregate_matrix
from rioxarray import open_rasterio
import numpy as np
import xarray
import matplotlib.pyplot as plt
from scipy.sparse import csr_matrix
from pandas import RangeIndex


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
        self.data = read_file(shapefile_path).set_index("location")
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
    cutout = Cutout(
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
            open_rasterio(raster_path, masked=True)
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

    The ambient temperatures are obtained directly from the atlite cutout,
    but the total effective irradiation `I_eff` is calculated as:

    .. math::
        I_{eff, horizontal} = I_{diff, horizontal} + c_{shading} I_{dir, horizontal} \\
        I_{eff, vertical} = I_{diff, vertical} + I_{ground, vertical} + \frac{c_{shading}}{\pi} I_{dir, vertical}

    where `I_diff` is the diffuse irradiation on the surface,
    `I_ground` is the ground reflected irradiation on the surface,
    `c_shading` is the external shading coefficient,
    and `I_dir` is the direct irradiation on the surface.
    The `I_dir,vertical` is calculated from atlite with vertical axis tracking,
    so it always matches the azimuth of the sun. Hence the division with
    pi to account for not all vertical surface area facing directly at the sun.

    Parameters
    ----------
    cutout : `atlite` cutout containing the weather data.
    external_shading_coefficient : `float`
        Coefficient for direct irradiation quantities to account for shading due to environment.

    Returns
    -------
    ambient_temperature_K : `xarray.DataSet`
        Ambient temperature timeseries [K].
    total_effective_irradiation_W_effm2 : `dict`
        Total effective irradiation timeseries [W/m2] for horizontal and vertical surfaces.
    """
    # Define the slope and azimuth angles for the horizontal and vertical surfaces.
    dirs = {  # Azimuths don't really matter due to slope and tracking respectively.
        "horizontal": get_orientation({"slope": 0.0, "azimuth": 0.0}),
        "vertical": get_orientation({"slope": 90.0, "azimuth": 0.0}),
    }

    # Get solar position and surface orientations
    solar_position = SolarPosition(cutout.data)
    surface_orientation = {
        dir: SurfaceOrientation(
            cutout.data,
            solar_position,
            orientation,
            "vertical",
        )
        for dir, orientation in dirs.items()
    }

    # Calculate total diffuse irradiation (diffuse + ground) and direct irradiation.
    diffuse_irradiation_W_m2 = {
        dir: TiltedIrradiation(
            cutout.data,
            solar_position,
            orientation,
            trigon_model="simple",
            clearsky_model="simple",
            tracking="vertical",
            irradiation="diffuse",
        )
        + TiltedIrradiation(
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
        dir: TiltedIrradiation(
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
    # Also note that the vertical direct irradiation is divided per `np.pi` to account for not all vertical surface area facing the right way!
    total_effective_irradiation_W_effm2 = {
        "horizontal": diffuse_irradiation_W_m2["horizontal"]
        + external_shading_coefficient * direct_irradiation_W_m2["horizontal"],
        "vertical": diffuse_irradiation_W_m2["vertical"]
        + external_shading_coefficient / np.pi * direct_irradiation_W_m2["vertical"],
    }

    # Return the weather quantities
    return cutout.data["temperature"], total_effective_irradiation_W_effm2


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

    Essentially, calculates the steady-state initial heating demand for
    the interior air of a building with the given parameters.
    Note that this demand is only part of the total demand,
    as it only accounts for heat gains and transfers applying
    directly to the indoor air of a building.

    The steady-state initial heating demand `\Phi_heating`
    for a given set-point on each time step is calculated as:

    .. math::
        \Phi_{heating} = \rho T' + \sigma (T' - T_a) - \phi_{int} - f_{conv} * f_{np} * g_{gl} * (A_{w,hor} * I_{eff,hor} + A_{w,ver} * I{eff,ver})

    where `\rho` is the self-discharge coefficient *(usually none)*,
    `\sigma` is the total heat transfer coefficient to the ambient air,
    `\phi_int` are the internal heat gains,
    `f_{conv}` is the assumed convective fraction of solar heat gains,
    `f_{np}` is the assumed non-perpendicularity correction factor for windows,
    `g_{gl}` is the assumed total normal solar energy transmittance of the glazing,
    `A_{w,hor}` is the horizontal window *(e.g. skylight)* surface area,
    `I_{eff,hor}` is the total effective solar irradiation on a horizontal surface,
    `A_{w,ver}` is the vertical window surface area,
    and `I_{eff,ver}` is the total effective solar irradiation on a vertical surface.
    Cooling demand is just the negative part of the steady-state
    heating demand.

    The above heating demand is only preliminary,
    as it only accounts for gains and losses of the indoor air node.
    Gains and losses through the envelope, as well as structural
    thermal mass dynamics are accounted for later in ArBuMo.jl.

    Parameters
    ----------
    set_point_K : xarray.DataArray
        The desired set point for the heating demand calculation.
    ambient_temperature_K : xarray.DataArray
        Ambient air temperature in Kelvin.
    total_effective_irradiation_W_effm2 : dict of xarray.DataArray
        Total effective irradiation in W per effective surface area,
        for horizontal and vertical surfaces separately.
    internal_heat_gains_W : xarray.DataArray
        Internal heat gains in W.
    self_discharge_coefficient_W_K : float
        Self-discharge coefficient in W/K.
    total_ambient_heat_transfer_coefficient_W_K : float
        Total heat transfer coefficient to ambient air in W/K.
    solar_heat_gain_convective_fraction : float
        Assumed convective fraction of the solar heat gains.
    window_non_perpendicularity_correction_factor : float
        Assumed window non-perpendicularity factor.
    total_normal_solar_energy_transmittance : float
        Assumed total normal solar energy transmittance of the windows.
    vertical_window_surface_area_m2 : float
        Vertical window surface area in m2.
    horizontal_window_surface_area_m2 : float
        Horizontal window (skylight) surface area in m2.

    Returns
    -------
    initial_heating_demand_W : xarray.DataArray
        The initial steady-state heating demand in W.
        Structural mass dynamics applied in ArBuMo.jl.
    """
    initial_heating_demand_W = (
        self_discharge_coefficient_W_K * set_point_K  # Self-discharge heat losses
        + total_ambient_heat_transfer_coefficient_W_K  # Ambient heat losses: This need to account for HRU bypass for cooling!
        * (set_point_K - ambient_temperature_K)  # Ambient heat losses.
        - internal_heat_gains_W  # NOTE! These need to include utilisable losses from the DHW tank!
        - solar_heat_gain_convective_fraction  # Solar gains are more complicated.
        * window_non_perpendicularity_correction_factor
        * total_normal_solar_energy_transmittance
        * (  # The effective irradiation already accounts for surface orientation.
            vertical_window_surface_area_m2
            * total_effective_irradiation_W_effm2["vertical"]
            + horizontal_window_surface_area_m2
            * total_effective_irradiation_W_effm2["horizontal"]
        )
    )
    return initial_heating_demand_W


def aggregate_xarray(xar, layout):
    """
    Aggregates xarray.DataArray data using a layout raster.

    Loosely based on what PyPSA/atlite `convert_and_aggregate` does.
    Unfortunately, that function couldn't be used as is due to the
    way it's programmed.

    Parameters
    ----------

    xar : xarray.DataArray
        The xarray data to be aggregated.
    layout : xarray.DataArray
        The layout used to spatially aggregate the `xar`.

    Results
    -------
    aggregate :
        The final aggregated data
    """
    # More or less copy-paste from PyPSA/atlite, no idea how this works.
    assert isinstance(layout, xarray.DataArray)
    layout = layout.reindex_like(xar).stack(spatial=["y", "x"])
    matrix = csr_matrix(layout.expand_dims("new"))
    index = RangeIndex(matrix.shape[0])
    return aggregate_matrix(xar, matrix=matrix, index=index)


def aggregate_demand_and_weather(
    shapefile_path,
    weather_start,
    weather_end,
    weights,
    external_shading_coefficient,
    heating_set_point_K,
    cooling_set_point_K,
    internal_heat_gains_W,
    self_discharge_coefficient_W_K,
    total_ambient_heat_transfer_coefficient_with_HRU_W_K,
    total_ambient_heat_transfer_coefficient_without_HRU_W_K,
    solar_heat_gain_convective_fraction,
    window_non_perpendicularity_correction_factor,
    total_normal_solar_energy_transmittance,
    vertical_window_surface_area_m2,
    horizontal_window_surface_area_m2,
    raster_path=None,
    resampling=5,
    filename="scope",
    save_layouts=True,
):
    """
    Calculates the aggregated initial heating/cooling demand and weather data.

    The master function of this module, handling the entire processing
    workflow from fetching ERA5 weather data to processing and aggregating
    the preliminary heating/cooling demands along with the weather for
    ArBuMo.jl.

    Parameters
    ----------
    shapefile : Shapefile
        The `Shapefile` used for defining the bounds of the cutout.
    weather_start : string
        The start of the cutout in `yyyy-mm-dd`, month and day can be skipped.
    weather_end: string
        The end of the cutout in `yyyy-mm-dd`, month and day can be skipped.
    weights : dict
        Weights for the `location_id` fields in the shapefile.
    external_shading_coefficient : float
        Coefficient for direct irradiation quantities to account for shading due to environment.
    heating_set_point_K : array
        The desired set points for the heating demand calculation as a time series.
    cooling_set_point_K : array
        The desired set points for the heating demand calculation as a time series.
    internal_heat_gains_W : array
        Internal heat gains in W, as a time series.
    self_discharge_coefficient_W_K : float
        Self-discharge coefficient in W/K.
    total_ambient_heat_transfer_coefficient_W_K : float
        Total heat transfer coefficient to ambient air in W/K.
    solar_heat_gain_convective_fraction : float
        Assumed convective fraction of the solar heat gains.
    window_non_perpendicularity_correction_factor : float
        Assumed window non-perpendicularity factor.
    total_normal_solar_energy_transmittance : float
        Assumed total normal solar energy transmittance of the windows.
    vertical_window_surface_area_m2 : float
        Vertical window surface area in m2.
    horizontal_window_surface_area_m2 : float
        Horizontal window (skylight) surface area in m2.

    Other parameters
    ----------------
    raster_path=None : Optional path to further raster data used for generating the layout.
    resampling=5 : Setting for resampling the data, `5=average` by default.
        See [`rasterio.enums.Resampling`](https://rasterio.readthedocs.io/en/stable/api/rasterio.enums.html#rasterio.enums.Resampling)
    filename = "scope" : str
        Name of the image file displaying the final layout.
    save_layouts = True : bool
        A flag to control whether layout images are to be saved.

    Returns
    -------
    """
    # Form the shapefile and prepare the cutout and layout
    shapefile = Shapefile(shapefile_path)
    cutout = prepare_cutout(shapefile, weather_start, weather_end)
    raster, layout = prepare_layout(shapefile, cutout, weights, raster_path, resampling)

    # Save layouts if desired.
    if save_layouts:
        for rst, name in [(raster, "raster"), (layout, "layout")]:
            f = plot_layout(shapefile, rst, dpi=600, title=filename + " " + name)
            f.savefig("figs/" + filename + "_" + name + ".png", dpi=600)
            plt.close(f)

    # Preprocess the weather
    ambient_temperature_K, total_effective_irradiation_W_effm2 = preprocess_weather(
        cutout, external_shading_coefficient
    )

    # Expand set points and internal heat gains to xarrays
    heating_set_point_K = expand_to_xarray(
        heating_set_point_K, ambient_temperature_K, "Heating set point", "K"
    )
    cooling_set_point_K = expand_to_xarray(
        cooling_set_point_K, ambient_temperature_K, "Heating set point", "K"
    )
    internal_heat_gains_W = expand_to_xarray(
        internal_heat_gains_W, ambient_temperature_K, "Internal heat gains", "W"
    )

    # Process initial heating demand, with HRU!
    heating_demand_W = process_initial_heating_demand(
        heating_set_point_K,
        ambient_temperature_K,
        total_effective_irradiation_W_effm2,
        internal_heat_gains_W,
        self_discharge_coefficient_W_K,
        total_ambient_heat_transfer_coefficient_with_HRU_W_K,
        solar_heat_gain_convective_fraction,
        window_non_perpendicularity_correction_factor,
        total_normal_solar_energy_transmittance,
        vertical_window_surface_area_m2,
        horizontal_window_surface_area_m2,
    )
    # Heating demand is indicated by the positive values of the steady-state solution,
    # cut out cooling demand as it likely has its own set point.
    heating_demand_W = xarray.where(heating_demand_W < 0.0, 0.0, heating_demand_W)

    # Process initial cooling demand, without HRU!
    cooling_demand_W = process_initial_heating_demand(
        cooling_set_point_K,
        ambient_temperature_K,
        total_effective_irradiation_W_effm2,
        internal_heat_gains_W,
        self_discharge_coefficient_W_K,
        total_ambient_heat_transfer_coefficient_without_HRU_W_K,
        solar_heat_gain_convective_fraction,
        window_non_perpendicularity_correction_factor,
        total_normal_solar_energy_transmittance,
        vertical_window_surface_area_m2,
        horizontal_window_surface_area_m2,
    )
    # Cooling demand is indicated by the negative values of the steady-state solution,
    # cut out heating demand as it likely has its own set point.
    cooling_demand_W = xarray.where(cooling_demand_W > 0.0, 0.0, -cooling_demand_W)

    # Aggregate and return the desired quantities.
    heating_demand_W = aggregate_xarray(heating_demand_W, layout)
    cooling_demand_W = aggregate_xarray(cooling_demand_W, layout)
    ambient_temperature_K = aggregate_xarray(ambient_temperature_K, layout)
    total_effective_irradiation_W_effm2["horizontal"] = aggregate_xarray(
        total_effective_irradiation_W_effm2["horizontal"], layout
    )
    total_effective_irradiation_W_effm2["vertical"] = aggregate_xarray(
        total_effective_irradiation_W_effm2["vertical"], layout
    )
    return (
        heating_demand_W,
        cooling_demand_W,
        ambient_temperature_K,
        total_effective_irradiation_W_effm2,
    )
