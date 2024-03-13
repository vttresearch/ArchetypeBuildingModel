# Input data processing for large-scale energy system modelling frameworks

This section aims to provide brief overviews of adapting the archetype building model
for use in large-scale energy system modelling frameworks like
[Backbone](https://cris.vtt.fi/en/publications/backbone)
and [SpineOpt](https://github.com/Spine-project/SpineOpt.jl).
In general, such frameworks are rarely designed with building modelling
in mind, and as such don't support the building-specific data stored within
[`ArchetypeBuilding`](@ref) and related `struct`s without further processing.
Perhaps the most important step is to understand the roles of
[`AbstractNode`](@ref) and [`AbstractProcess`](@ref) as opposed to their
detailed counterparts [`BuildingNodeData`](@ref) and [`BuildingProcessData`](@ref).

The following [Processing thermal nodes into `AbstractNode`s](@ref)
and [Processing HVAC equipment into `AbstractProcess`es](@ref) sections
explain the aggregation and abstraction of building level data for energy-system-scale,
while the  [Backbone input data processing](@ref) and [SpineOpt input data processing](@ref)
sections focus on the processing workflows of their respective models models.


### Processing thermal nodes into `AbstractNode`s

[`BuildingNodeData`](@ref) has been designed to be as human-readable
and understandable as possible from a building-domain point-of-view,
including a plethora of different parameters related to archetype definitions,
assumptions, sizing, etc.
However, most of these parameters don't have counterparts in large-scale energy
system models, and need to be aggregated into something meaningful.
This is done via the [`ArBuMo.process_abstract_node`](@ref) function,
aggregating all the separate categories of thermal mass, heat transfer,
as well as ambient condition and external load related parameters
into the bare essentials required for depicting the node in large-scale
energy system modelling frameworks, represented by an [`AbstractNode`](@ref).

!!! note 
    Since ambient-temperature-related interactions are rarely directly supported by large-scale energy system modelling frameworks, they are re-cast into *self-discharge* and *external load* components as explained in the [`ArBuMo.process_abstract_node`](@ref) docstring.

As this abstraction reduces the amounts of different terms in the equations,
it also happens to simplify
[Solving the baseline heating demand and HVAC equipment consumption](@ref),
which is why the [`AbstractNode`](@ref)s are used for the calculations in the
[`ArBuMo.solve_heating_demand`](@ref) function.


## Backbone input data processing

This section aims to provide an overview of the processing done for
producing [Backbone](https://cris.vtt.fi/en/publications/backbone) input data,
explaining the logic of and functions withing `src/create_backbone_input.jl`.

The [`BackboneInput`](@ref) `struct` contains the relevant input data structure
for representing the produced archetype building models
within [Backbone](https://cris.vtt.fi/en/publications/backbone).
The similarly named constructor takes as input a dictionary containing the
defined [building\_archetype](@ref) linked to its processed
[`ArchetypeBuildingResults`](@ref),
loops over the archetypes, and adds them to the [`BackboneInput`](@ref)
one by one.
The actual processing is handled by the
[`ArBuMo.add_archetype_to_input!`](@ref) function,
which performs a lot of rather complicated manipulations to adapt the archetype
building data for [Backbone](https://cris.vtt.fi/en/publications/backbone).
For people familar with [Backbone](https://cris.vtt.fi/en/publications/backbone)
model structure, the key points are:

1. Each [building\_archetype](@ref) is mapped into a `grid`.
2. Each [`AbstractNode`](@ref) is mapped into a `node` in the corresponding archetype `grid`.
2. Each [`AbstractProcess`](@ref) in each [building\_archetype](@ref) is mapped into a unique `unit`.
4. System link nodes defined by [building\_archetype\_\_system\_link\_node](@ref) relationships and the associated [node\_name](@ref) and [grid\_name](@ref) are created to serve as connection points to potential energy system datasets.

Since [ArBuMo.jl](@ref) is based on
[SpineInterface.jl](https://github.com/Spine-project/SpineInterface.jl),
the produced [Backbone](https://cris.vtt.fi/en/publications/backbone) input
is saved in its Spine Datastore format, and requires the use of
[Spine Toolbox](https://github.com/Spine-project/Spine-Toolbox) and the
associated tools contained within the [Backbone](https://cris.vtt.fi/en/publications/backbone)
repository in order to produce the `inputData.gdx` for running the model.
Additionally, the `backbone_utils/export_auxiliary_building_data.json` exporter
specification can be used for exporting useful `.gdx` data not directly used
by [Backbone](https://cris.vtt.fi/en/publications/backbone).


## SpineOpt input data processing

This section aims to provide and overview of the processing done for
producing [SpineOpt](https://github.com/Spine-project/SpineOpt.jl)
input data, explaining the logic and functions within
`src/create_spineopt_input.jl`.

The [`SpineOptInput`](@ref) `struct` contains the relevant input data structure
for representing the produced archetype building models within
[SpineOpt](https://github.com/Spine-project/SpineOpt.jl).
The similarly named constructor takes as input a dictionary containing the
defined [building\_archetype](@ref) linked to its processed
[`ArchetypeBuildingResults`](@ref),
loops over the archetypes, and adds them to the [`SpineOptInput`](@ref)
one by one.
The actual processing is handled by the
[`ArBuMo.add_archetype_to_input!`](@ref) function,
which essentially maps the [`AbstractNode`](@ref) and [`AbstractProcess`](@ref)
parameters to their [SpineOpt](https://github.com/Spine-project/SpineOpt.jl) 
counterparts. For people familiar with
[SpineOpt](https://github.com/Spine-project/SpineOpt.jl),
the key points are:

1. Each [`AbstractNode`](@ref) in each [building\_archetype](@ref) is mapped into a unique `node`.
2. Each [`AbstractProcess`](@ref) in each [building\_archetype](@ref) is mapped into a unique `unit`.
3. System link nodes defined by [building\_archetype\_\_system\_link\_node](@ref) relationships and the associated [node\_name](@ref) are created to serve as connection points to potential energy system datasets.

Since both [ArBuMo.jl](@ref) and
[SpineOpt](https://github.com/Spine-project/SpineOpt.jl) are based on
[SpineInterface.jl](https://github.com/Spine-project/SpineInterface.jl),
the produced SpineOpt input data is immediately compatible.
However, note that the produced SpineOpt input data only contains the
data describing the modelled building stock, without all the necessary
definitions to run the model in any meaningful way.
