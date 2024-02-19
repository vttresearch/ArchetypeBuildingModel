# ArBuMo.jl

[![license](https://img.shields.io/badge/license-MIT-brightgreen)](https://mit-license.org/)

A Julia module for aggregating building stock data into desired synthetic average archetype building
lumped-capacitance thermal models for large-scale energy system modelling purposes.

The goal of this module is to aggregate building stock data into arbitrary
simplified lumped-capacitance thermal models representing the flexible
heating/cooling demand of the underlying building stock for large-scale
energy system models like [Backbone](https://cris.vtt.fi/en/publications/backbone)
or [SpineOpt](https://github.com/Spine-project/SpineOpt.jl).
Requires building stock data processed by e.g. [FinnishBuildingStockData.jl](https://github.com/vttresearch/FinnishBuildingStockData) or similar,
as well as archetype building definitions according to the template provided in `archetype_definitions.json` as input.
For automatic weather data processing, this module uses the Python sub-module `ArBuWe`.

ArBuMo.jl is a fork of [ArchetypeBuildingModel.jl](https://github.com/vttresearch/ArchetypeBuildingModel) for further developments,
as implementing investment modelling options requires reformulating the building stock heating/cooling demand modelling
as well as the lumped-capacitance thermal models.

>[!IMPORTANT]
>This repository is currently undergoing major changes, and a lot of the information in this README
>and in the documentation in particular is likely to be outdated.
>Hopefully I'll be done with most of the revisions sometime in late Spring 2024,
>and can focus on making sure all presented information is up to date.


## Key contents

1. `archetype_definitions.json`, a template for the Spine data structure required for defining the archetype buildings.
2. `process_archetype_buildings.jl`, the main program file for archetype building processing using Spine Toolbox.
3. `process_archetype_buildings.json`, the [Spine Toolbox](https://github.com/Spine-project/Spine-Toolbox) tool specification for the main program above.
4. `testscript.jl`, a testing script/tutorial for using this module.
5. `testscript.ipynb`, a testing notebook/tutorial for using the `archetypebuildingweather` sub-module.
6. `src/ArBuMo.jl`, the main Julia module file.
7. `src/archetypebuildingweather/` contains the Python sub-module for automatic weather data processing.
8. `data/` contains data related to the automatic weather data processing.
9. `examples/` contains example definitions for `building_fabrics`, `building_systems`, and `building_weather`.
10. `figs/` contains automatically output diagnostic figures about the automatic weather data aggregation. 


## ArBuWe.py

>[!IMPORTANT]
>The weather data processing will need to be heavily revised in order to enable investment modelling options.
>Everything explained below is subject to major changes.

A Python sub-module for processing weather data.

This package contains the necessary functions for aggregating ERA5 weather data
using [`PyPSA/atlite`](https://github.com/PyPSA/atlite) for large geographical
areas described using a shapefile, a set of weights connected to the
shapefile, and an optional raster data file also used for weighting.
Essentially, the operating workflow of the module goes as follows:

1. Load a *shapefile* defining the scope of weather data for `atlite`, and fetch the data.
2. Assign the given *weights* to the *shapefile*, and include the optional *raster* data if defined.
3. Downsample and normalize the *weight raster* generated above to match ERA5 resolution.
4. Estimate the aggregated steady-state raw heating/cooling demand for the indoor air node for each archetype building and heating/cooling set-point.
5. Calculate aggregated weather parameters required for `ArBuMo.jl` for estimating the final total heating/cooling demand when accounting for thermal mass dynamics.


## Installation

Unfortunately, this module has quite complicated dependencies due to
[`SpineInterface.jl`](https://github.com/Spine-project/SpineInterface.jl) and
[`PYPSA/atlite`](https://github.com/PyPSA/atlite) needing to operate using a
single python environment accessible via [`PyCall`](https://github.com/JuliaPy/PyCall.jl).
Thus, it needs quite a specific python environment to be set up and linked with
[`SpineInterface.jl`](https://github.com/Spine-project/SpineInterface.jl).
Here's a rough sketch of how I got it working:

1. Create a conda environment with python 3.8, e.g. `conda create -n toolboxatlite python=3.8` and activate it.
2. Install [Spine Toolbox](https://github.com/Spine-project/Spine-Toolbox) from a local source via `pip install -r requirements.txt` in the Toolbox root folder.
3. Install `PYPSA/atlite` via `conda install -c conda-forge atlite`
4. Install the `ArchetypeBuildingWeather.py` sub-module via `pip install -e .` in the `src/` folder.
5. Install JupyterLab via `conda install -c conda-forge jupyterlab`.


## Usage

This module is intended to be used as a part of a [Spine Toolbox](https://github.com/Spine-project/Spine-Toolbox) workflow,
with the following rough steps:

1. Create a Spine Datastore with building stock data processed by e.g. [FinnishBuildingStockData.jl](https://github.com/vttresearch/FinnishBuildingStockData).
2. Import `archetype_definitions.json` on top of the processed building stock data.
3. Define the desired `building_archetype` objects, and connect them to the relevant `building_fabrics`, `building_loads`, `building_scope`, and `building_system` objects. *(`building_weather` is optional, and will be automatically generated if missing.)*
4. Use the `process_archetype_buildings.json` tool to process the data and definitions into the desired energy system model input.

In case you need to familiarize yourself with the inner workings of this module,
see `testscript.jl` for examples on how to use the `ArBuMo.jl`,
and `testscript.ipynb` for examples on how to use the `ArBuWe.py` sub-module.
Note that it is recommended to store the relevant geographical information system (GIS)
data used via the module within the `data/` folder in the repository for the moment.


## Documentation

>[!NOTE]
>Online documentation is temporarily offline, as there are major reformulations pending to enable investment modelling options.
>The documentation is likely not up to date as a result!

In order to build and read the documentation locally,
start a Julia REPL from the root folder of this module and perform the following steps:

1. Activate the `documentation` environment from the Julia Package manager
```julia
julia> ]
(ArBuMo) pkg> Activate documentation
(documentation) pkg> ]
julia>
```

2. Run the `documentation/make.jl` script to build the documentation.
```julia
julia> include("documentation/make.jl")
```

3. Open the newly built `documentation/build/index.html` to start browsing the documentation.


## License

MIT, see `LICENSE` for more information.


## How to cite

TBD


## Acknowledgements

<center>
<table width=500px frame="none">
<tr>
<td valign="middle" width=100px>
<img src=https://www.aka.fi/globalassets/vanhat/y_kuvat/aka_logo_en.svg alt="AKA emblem" width=100%></td>
<td valign="middle">
This module was built for the Academy of Finland project "Integration of building flexibility into future energy systems (FlexiB)" under grant agreement No 332421.
</td>
</table>
</center>

<center>
<table width=500px frame="none">
<tr>
<td valign="middle" width=100px>
<img src=https://european-union.europa.eu/themes/contrib/oe_theme/dist/eu/images/logo/standard-version/positive/logo-eu--en.svg alt="EU emblem" width=100%></td>
<td valign="middle">
This work has been partially supported by EU project Mopo (2023-2026), which has received funding from European Climate, Infrastructure and Environment Executive Agency under the European Union’s HORIZON Research and Innovation Actions under grant agreement N°101095998.
</td>
</table>
</center>
