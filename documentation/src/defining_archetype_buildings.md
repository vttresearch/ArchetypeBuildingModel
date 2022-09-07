# Defining archetype buildings

This section aims to go over the key components of the
[Input data reference](@ref) and explain how these components
are used to define the archetype building lumped-capacitance thermal models.
The explanations remain at a pretty high level, but links to the
appropriate sections of the [Input data reference](@ref)
and [Library](@ref) sections are provided for readers interested in the
technical details.


## The `building_archetype` definition

**The single most important object class in the definitions**, as each
[building\_archetype](@ref) object essentially represents a single synthetic
average archetype building lumped-capacitance thermal model to be created.
In the end, the rest of the object classes used for defining different aspects
of the archetype buildings are each connected to a [building\_archetype](@ref)
using the appropriate relationship classes:
- [building\_archetype\_\_building\_scope](@ref): Defines one [building\_scope](@ref) for the [building\_archetype](@ref). Determines the geographical and statistical scope, e.g. the number, gross-floor area, and structural properties of the archetype based on aggregating the underlying building stock statistics.
- [building\_archetype\_\_building\_fabrics](@ref): Defines one [building\_fabrics](@ref) for the [building\_archetype](@ref). Determines how the different types of structures and the indoor air and furniture are aggregated into a lumped-capacitance thermal model depicting the archetype.
- [building\_archetype\_\_building\_systems](@ref): Defines one [building\_systems](@ref) for the [building\_archetype](@ref). Determines the HVAC system of the archetype.
- [building\_archetype\_\_building\_loads](@ref): Defines one [building\_loads](@ref) for the [building\_archetype](@ref). Determines the internal loads and domestic hot water demand for the archetype.
- [building\_archetype\_\_building\_weather](@ref): Defines which [building\_weather](@ref) to use for the [building\_archetype](@ref). If left undefined, the [ArchetypeBuildingWeather.py](@ref) sub-module will be used to try and automatically fetch and process weather data based on the [building\_scope](@ref) and other parameters.
- [building\_archetype\_\_system\_link\_node](@ref): Allows customizable definition of nodes intended to be used as the connection points between the lumped-capacitance thermal models and the overarching large-scale energy system models.

The [building\_archetype](@ref) objects also houses most of the parameters
defining assumptions regarding the archetype and how it's modelled.
The only mandatory parameter is the [weather\_year](@ref), which is required
for automatic weather data processing using [ArchetypeBuildingWeather.py](@ref),
while the rest of the parameters have set default values that kick in if not
specified by the user.
For example, the assumed shape of the building envelope can be tweaked using
e.g. [building\_frame\_depth\_m](@ref), [number\_of\_storeys](@ref),
[room\_height\_m](@ref), and [window\_area\_to\_external\_wall\_ratio\_m2\_m2](@ref).
Similarly, some assumptions related to the modelling can be tweaked here,
e.g. [effective\_thermal\_capacity\_of\_interior\_air\_and\_furniture\_J\_m2K](@ref),
[external\_shading\_coefficient](@ref), or [internal\_heat\_gain\_convective\_fraction](@ref).
See the [Object parameters](@ref) section for a comprehensive list of all the
possible parameters for the [building\_archetype](@ref).

Ultimately, each [building\_archetype](@ref) is processed into an [`ArchetypeBuilding`](@ref)
by the main program described by the [Overview of the workflow](@ref) section.


## The `building_scope` definition

Arguably the second-most important object class in the definiton,
as each [building\_scope](@ref) object describes how to aggregate the
underlying building stock statistics to calculate the properties of a
synthetic average archetype building.
Essentially, the [building\_scope](@ref) object tells which [building\_stock](@ref)
data to use, and which [building\_type](@ref)s, [heat\_source](@ref)s,
and [location\_id](@ref)s to include, and from what period of time.
Similar to how [The `building_archetype` definition](@ref) works,
the [building\_scope](@ref) is also defined largely by connecting the desired
pieces to it using the appropriate relationship classes:
- [building\_archetype\_\_building\_scope](@ref): Defines one [building\_scope](@ref) for the [building\_archetype](@ref). Determines the geographical and statistical scope, e.g. the number, gross-floor area, and structural properties of the archetype based on aggregating the underlying building stock statistics.
- [building\_scope\_\_building\_stock](@ref): Defines which [building\_stock](@ref) data to use, and with what weights. Allows e.g. pseudo-interpolating data for years between recorded [building\_stock](@ref) years.
- [building\_scope\_\_building\_type](@ref): Defines which [building\_type](@ref)s to include and are represented by the [building\_scope](@ref), with user-defined weights if desired.
- [building\_scope\_\_heat\_source](@ref): Defines which [heat\_source](@ref)s to include and are represented by the [building\_scope](@ref), with user-defined weights if desired.
- [building\_scope\_\_location\_id](@ref): Defines which [location\_id](@ref)s to include and are represented by the [building\_scope](@ref), with user-defined weights if desired.

Unlike the rest of the geographical and statistical scope definitions,
the temporal aggregation is handled via parameters instead of relationships.
This is done to facilitate selection of the desired time periods without
having to look into the underlying data to know of the exact time periods
used in the statistical data.
Similarly, the years for the ventilation and structural properties rarely
match the time periods for the building stock statistics, and can be sampled
as easily using parameters instead.
For this purpose, the [building\_scope](@ref) has the two mandatory parameters:
- [scope\_period\_start\_year](@ref): The first year included in the [building\_scope](@ref).
- [scope\_period\_end\_year](@ref): The last year included in the [building\_scope](@ref).

Ultimately, each [building\_scope](@ref) is processed into a [`ScopeData`](@ref)
by the main program described by the [Overview of the workflow](@ref) section.


## The `building_fabrics` definition

The [building\_fabrics](@ref) objects are used to define how the [building\_archetype](@ref)
envelope, interior structures, and interior air and furniture are modelled.
Essentially, it defines the temperature nodes of the structural components of
the building for the lumped-capacitance thermal model of the archetype.
Different [building\_fabrics](@ref) are relatively lightweight to define,
but they rely heavily on [The `building_node` definition](@ref)s,
which are discussed in the following subsection.
Regardless, [building\_fabrics](@ref) are simply collections of
[building\_node](@ref) using the following relationship classes:
- [building\_archetype\_\_building\_fabrics](@ref): Defines one [building\_fabrics](@ref) for the [building\_archetype](@ref). Determines how the different types of structures and the indoor air and furniture are aggregated into a lumped-capacitance thermal model depicting the archetype.
- [building\_fabrics\_\_building\_node](@ref): Defines the [building\_node](@ref)s constituting this [building\_fabrics](@ref).

The [building\_fabrics](@ref) definitions help form the [`EnvelopeData`](@ref)
during the data processing by the main program, described by the
[Overview of the workflow](@ref) section.


### The `building_node` definition

The [building\_node](@ref) objects are used to define the temperature nodes
for the lumped-capacitance thermal models, whether they be a part of the
[building\_fabrics](@ref), e.g. building envelope structures,
or the [building\_systems](@ref), e.g. a domestic hot water storage tank.
Due to their rather abstract nature, [building\_node](@ref)s have quite a few
different ways they can be defined, depending on the desired purpose.
Regardless, the principles remain similar to above, and most of the
[building\_node](@ref) definitions are handled via the following relationship
classes:
- [building\_archetype\_\_system\_link\_node](@ref): Designates whether the [building\_node](@ref) is a special *system link node*, meaning that it is used as the connection point between the generated lumped-capacitance thermal model and the overarching large-scale energy models. Note that only the [@system\_link\_node\_1](@ref) etc. should ever be used for this relationship class.
- [building\_fabrics\_\_building\_node](@ref): Defines the [building\_node](@ref)s constituting the [building\_fabrics](@ref).
- [building\_node\_\_building\_node](@ref): Allows user-defined custom heat transfer coefficients between [building\_node](@ref)s. Should only be necessary for more advanced use-cases.
- [building\_node\_\_structure\_type](@ref): Defines which [structure\_type](@ref)s *(if any)* are included on this [building\_node](@ref), affecting its effective thermal mass and heat transfer coefficients between the interior and exterior of the archetype. User-defined [structure\_type\_weight](@ref)s can be provided to tweak the model further, but shouldn't be necessary for basic use-cases.
- [building\_process\_\_direction\_\_building\_node](@ref): Defines how [building\_process](@ref)es interact with the [building\_node](@ref)s, and houses their maximum power flow parameters.
- [building\_systems\_\_building\_node](@ref): Defines the [building\_node](@ref)s included in the [building\_systems](@ref).

The [building\_node](@ref) object class also contains a few important parameters.
A comperehensive list can be found in the [Object parameters](@ref) section,
but the most important ones are:
- [domestic\_hot\_water\_demand\_weight](@ref): Defines the share of the total domestic hot water demand attributed to this [building\_node](@ref). Typically set to either zero or one, but e.g. if the DHW tank is modelled using several nodes, tweaking the value might be necessary.
- [interior\_air\_and\_furniture\_weight](@ref): Defines the share of the total interior air and furniture attributed to this [building\_node](@ref), affecting its effective thermal mass, heat transfer coefficients due to ventilation/infiltration and windows, etc. Typically set to either zero or one, but if multi-zone modelling is attempted *(not recommended)*, tweaking the value might be necessary.
- [maximum\_permitted\_temperature\_K](@ref) and [minimum\_permitted\_temperature\_K](@ref) are used to set the temperature limits for the node.

Ultimately, each [building\_node](@ref) is processed into a [`BuildingNodeData`](@ref)
*(and an [`AbstractNode`](@ref) for energy system model export)*
by the main program described by the [Overview of the workflow](@ref) section.


## The `building_systems` definition

The [building\_systems](@ref) objects are used to define different heating/cooling
systems for the archetype buildings. They are similar to [building\_fabrics](@ref)
in the sense that they rely heavily on more detailed definitions for
[building\_node](@ref)s and [building\_process](@ref)s via these
relationship classes:
- [building\_archetype\_\_building\_systems](@ref): Defines one [building\_systems](@ref) for the [building\_archetype](@ref). Determines the HVAC system of the archetype.
- [building\_systems\_\_building\_node](@ref): Defines [building\_node](@ref)s parts of this [building\_systems](@ref), e.g. a domestic hot water storage tank.
- [building\_systems\_\_building\_process](@ref): Defines [building\_process](@ref)s parts of this [building\_systems](@ref), e.g. a direct electric resistance heater or a heat pump.
[The `building_node` definition](@ref) is already explained in the sections above,
but [building\_process](@ref)es are extremely important for [building\_systems](@ref)s.


### The `building_process` definition

The [building\_process](@ref) objects are used to define energy transfer/transformation
processes in the [building\_systems](@ref), e.g. direct electric resistance heaters
or heat pumps. The important relationship classes for the definitions are:
- [building\_process\_\_direction\_\_building\_node](@ref): Defines how [building\_process](@ref)es interact with the [building\_node](@ref)s, and houses their maximum power flow parameters.
- [building\_systems\_\_building\_process](@ref): Defines [building\_process](@ref)s parts of this [building\_systems](@ref). Also houses the maximum power flow parameters for the process.

The [building\_process](@ref) object class also contains a few important parameters,
mostly due to the potential complexity of different heat pump systems.
A comperehensive list can be found in the [Object parameters](@ref) section,
but the most important ones are:
- [coefficient\_of\_performance\_base](@ref): The base COP of the process.
- [coefficient\_of\_performance\_mode](@ref): Wherther the [building\_process] is used for `:heating` or `:cooling`.
See the [`ArchetypeBuildingModel.calculate_cop`](@ref) function for more details on how weather-dependent
COPs are modelled.

Ultimately, each [building\_process](@ref) is processed into a [`BuildingProcessData`](@ref)
*([`AbstractProcess`](@ref) for the energy sysem model export)*
by the main program described by the [Overview of the workflow](@ref) section.


## The `building_loads` definition

The [building\_loads](@ref) objects are used to define the internal heat loads
and domestic hot water demand for a [building\_archetype](@ref).
The only relationship class necessary is the: 
- [`building\_archetype\_\_building\_loads](@ref): Defines one [building\_loads](@ref) for the [building\_archetype](@ref).

For [building\_loads](@ref), all the important information is contained in the
[Object parameters](@ref):
- [domestic\_hot\_water\_demand\_base\_W](@ref): Fixed DHW demand data.
- [domestic\_hot\_water\_demand\_gfa\_scaling\_W\_m2](@ref): Gross-floor area scaling DHW demand data.
- [internal\_heat\_loads\_base\_W](@ref): Fixed internal heat gains data.
- [internal\_heat\_loads\_gfa\_scaling\_W\_m2](@ref): Gross-floor area scaling internal heat gains data.

Ultimately, each [building\_loads](@ref) is processed into a [`LoadsData`](@ref)
by the main program described by the [Overview of the workflow](@ref) section.


## The `building_weather` definition

The [building\_weather](@ref) objects are used for defining weather data for
a [building\_archetype](@ref). Of all the definitions described in this section,
it is the only non-required one. In case the relationship:
- [building\_archetype\_\_building\_weather](@ref): Defines which [building\_weather](@ref) to use for the [building\_archetype](@ref).
is left undefined, the [ArchetypeBuildingWeather.py](@ref) sub-module will be used
in an attempt to automatically fetch and process the weather data necessary for
creating the archetype lumped-capacitance thermal model.

In case the user wants to specify the weather data exactly, the following
[Object parameters](@ref) are required:
- [ambient\_temperature\_K](@ref): The ambient temperature in [K].
- [diffuse\_solar\_irradiation\_W\_m2](@ref): The diffuse horisontal irradiation in [W/m2].
- [direct\_solar\_irradiation\_W\_m2](@ref): A `Map` specifying the direct irradiation for vertical surfaces facing in all four cardinal directions in [W/m2].