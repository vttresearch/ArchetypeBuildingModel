# ArchetypeBuildingModel.jl

A Julia module for aggregating building stock data into desired
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
[Overview of the model workflow](@ref) section, and produces
lumped-capacitance thermal models of the desired synthetic average archetype
buildings depicting the heating/cooling demand and HVAC energy consumption
of a building stock,
explained in the [Archetype building modelling](@ref) section.
While weather data can be provided by the user if desired,
the [ArchetypeBuildingWeather.py](@ref) python sub-module provides automatic
fetching and processing of the necessary weather data based on the provided
archetype building definitions.
A baseline heating/cooling demand and HVAC consumption is also calculated
using very simple rule-based control aiming to keep the node temperatures
within permitted limits.
The key outputs from this module, however, are the readily made
[Backbone](https://cris.vtt.fi/en/publications/backbone) or
[SpineOpt](https://github.com/Spine-project/SpineOpt.jl) input datasets
that can be plugged into their respective energy system models for
depicting the flexible heating/cooling demand of the depicted building stock.

This documentation is organized as follows:
The [Defining archetype buildings](@ref) section explains how the archetype
buildings are defined, meaning the key components in the
[Input data reference](@ref), and how to use them.
The [Overview of the model workflow](@ref) section goes through the
`process_archetype_buildings.jl` main program file, explaining what is
actually being done when aggregating the building stock data into the
desired synthetic average archetype buildings.
The [Archetype building modelling](@ref) section explains the lumped-capacitance
thermal modelling approach used by this module in more detail,
while the [ArchetypeBuildingWeather.py](@ref) section briefly explains
the logic and workings of the python sub-module handling the automatic
weather data processing.
Finally, the [Input data reference](@ref) and the [Library](@ref)
sections provide comprehensive documentation of the definition/input
data format and the modelling code respectively.


## Contents

```@contents
```