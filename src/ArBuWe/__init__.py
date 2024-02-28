# __init__.py

# Initializes the `ArBuWe` module.

import logging

logging.basicConfig(level=logging.INFO)

from .processing import Shapefile
from .processing import prepare_cutout
from .processing import prepare_layout
from .processing import plot_layout
from .processing import preprocess_weather
from .processing import expand_to_xarray
from .processing import process_initial_heating_demand
from .processing import aggregate_xarray
from .processing import aggregate_demand_and_weather
