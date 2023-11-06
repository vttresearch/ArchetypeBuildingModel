# ArchetypeBuildingModel.jl

A [SpineInterface.jl](https://github.com/Spine-project/SpineInterface.jl)-based
Julia module for aggregating building stock data into desired
archetype building lumped-capacitance thermal models.

The goal of this module is to provide an easy way to define and aggregate
building stock statistical data into arbitrary sets of synthetic average
archetype building lumped-capacitance thermal models, depicting the flexible
heating/cooling demand of a building stock.
These lumped-capacitance thermal models are created primarily for
seamless integration into large-scale energy system models like
[Backbone](https://cris.vtt.fi/en/publications/backbone) or
[SpineOpt](https://github.com/Spine-project/SpineOpt.jl),
in order to depict flexible heating/cooling demand of significant portions
of the building stock.

Essentially, this module takes input data and archetype building definitions
in the format detailed in the [Input data reference](@ref)
section as input, processes them according to the workflow detailed in the
[Overview of the workflow](@ref) section, and produces
lumped-capacitance thermal models of the desired synthetic average archetype
buildings depicting the heating/cooling demand and HVAC energy consumption
of a building stock, explained in the [Archetype building modelling](@ref) section.
While weather data can be provided by the user if desired,
the [ArchetypeBuildingWeather.py](@ref) python sub-module provides automatic
fetching and processing of the necessary weather data based on the provided
archetype building definitions.
[Solving the baseline heating demand and HVAC equipment consumption](@ref)
is also calculated using very simple rule-based control keeping the node
temperatures within permitted limits.

The key outputs from this module, however, are the readily made
[Backbone](https://cris.vtt.fi/en/publications/backbone) or
[SpineOpt](https://github.com/Spine-project/SpineOpt.jl) input datasets
that can be plugged into their respective energy system models for
depicting the flexible heating/cooling demand of the depicted building stock.
See the figure below for a high-level illustration of the workflow.

This documentation is organized as follows:
The [Defining archetype buildings](@ref) section explains how the archetype
buildings are defined, meaning the key components in the
[Input data reference](@ref), and how to use them.
The [Overview of the workflow](@ref) section goes through the
`process_archetype_buildings.jl` main program file, explaining what is
actually being done when aggregating the building stock data into the
desired synthetic average archetype buildings.
The [Archetype building modelling](@ref) section explains the lumped-capacitance
thermal modelling approach used by this module in more detail,
while the [ArchetypeBuildingWeather.py](@ref) section briefly explains
the logic and workings of the python sub-module handling the automatic
weather data processing.
The [Input data processing for large-scale energy system modelling frameworks](@ref)
section provides an overview of how the data is further processed to be compatible
with [Backbone](https://cris.vtt.fi/en/publications/backbone)
and [SpineOpt](https://github.com/Spine-project/SpineOpt.jl).
Finally, the [Input data reference](@ref) and the [Library](@ref)
sections provide comprehensive documentation of the definition/input
data format and the modelling code respectively.

![ABMWorkflow](ABMFlow.png)


## Related works

For an up-to-date list of works using [ArchetypeBuildingModel.jl](@ref),
please refer to the
[VTT Research Information System entry on ArchetypeBuildingModel.jl](https://cris.vtt.fi/en/publications/archetypebuildingmodeljl-a-julia-module-for-aggregating-building-).
Regardless, here are a couple of author-curated highlights:

1. [Sensitivity analysis and comparison against dedicated white-box building simulation software (Preprint)](https://zenodo.org/doi/10.5281/zenodo.7623739)
2. [Building-level optimal energy management (Preprint)](https://zenodo.org/doi/10.5281/zenodo.7767363)
3. [Stochastic district-level energy management (Accepted and presented at the BS2023 conference, yet to be published)](https://cris.vtt.fi/en/publications/stochastic-model-predictive-control-of-district-scale-building-en)
4. [Demo for generating EU-level flexible building stock heating/cooling models for SpineOpt](https://zenodo.org/doi/10.5281/zenodo.8238141)
5. [Analysing the impacts of flexible residential electric heating for Finland (MSc thesis)](https://cris.vtt.fi/en/publications/heating-demand-response-in-detached-houses-comparing-cost-savings)


## Contents

```@contents
```