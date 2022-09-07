# Overview of the workflow

The following sections aim to go over the main program file `process_archetype_buildings.jl`
and explain the high-level workflow of the `ArchetypeBuildingModel.jl`.
Links are provided to the detailed explanations of the different data structures
and functions in the [Library](@ref) section for readers interested in the
technical details.


## Command line arguments

The `process_archetype_buildings.jl` main program file has been primarily
designed to be run via [Spine Toolbox](https://github.com/Spine-project/Spine-Toolbox),
but be run directly from the command line as well if necessary.
Regardless, the main program is controlled using the following command line arguments:
1. The url to a *Spine Datastore* containing the required input data and archetype building definitions.
Furthermore, the following optional keyword arguments can be provided:
- `-spineopt <url>`, the url to a *Spine Datastore* where the produced *SpineOpt* input data should be written, if any.
- `-backbone <url>`, the url to a *Spine Datastore* where the produced *Backbone* input data should be written, if any.
- `-import_weather <false>`, controls whether auto-generated `building_weather` are imported into the input *Datastore*. Set to `false` by default. 
- `-save_layouts <false>`, controls whether auto-generated `building_weather` layouts are saved as images. Set to `false` by default.

Essentially, the main purpose of the `-import_weather` keyword is to avoid
reprocessing the weather data over and over again for cases where multiple
simulations are run without changing the weather data.
Similarly, the `-save_layouts` keyword exists primarily for debugging purposes,
allowing visual inspection of the weather data weighting rasters.


## Input database tests

Next, the main program opens the input *Datastore*, and performs a series
of tests on the input data and archetype building definitions to check that
they make sense. If not, the main program will halt and display the test
results and error messages for the user, in order to help them deduce what is
wrong with the input data or definitions.

The input *Datastore* tests are handled by the [`run_object_class_tests`](@ref),
[`run_parameter_tests`](@ref), and [`run_structure_type_tests`](@ref) functions.


## Process `ScopeData` structs

As explained by [The `building_scope` definition](@ref) section,
the [building\_scope](@ref) defines the geographical and statistical scope
for a [building\_archetype](@ref). Before we can begin creating
lumped-capacitance thermal models for the archetype buildings, the main program
first needs to know the aggregated average basic properties of the archetype.
Thus, the next step is to process and create the [`ScopeData`](@ref) structs
for all the [building\_scope](@ref)s attached to a [building\_archetype](@ref)
via a [building\_archetype\_\_building\_scope](@ref) relationship
*(scopes not attached to any archetype are not processed)*.

The final [`ScopeData`](@ref) structs are stored into a
`scope_data_dictionary`, which can be examined through the Julia REPL
after the main program has finished.


## Process `WeatherData` structs

After the main program is done with [Process `ScopeData` structs](@ref),
the next step is to process the [building\_weather](@ref) definitions into
[`WeatherData`](@ref) before forming the archetype building lumped-capacitance
thermal models. However, since [The `building_weather` definition](@ref)
isn't mandatory, the main program first checks which [building\_archetype](@ref)
definitions lack a [building\_archetype\_\_building\_weather](@ref) definition,
and tries to fetch it automatically using the [ArchetypeBuildingWeather.py](@ref)
sub-module. Essentially, the automatic weather generation is handled via the
[`create_building_weather`](@ref) function based on
[The `building_archetype` definition](@ref) and the processed [`ScopeData`](@ref).

If the `-import_weather true` argument has been given, the main program will
attempt to import the automatically generated [building\_weather](@ref) object
into the input and definition *Datastore*, as well as link it to the appropriate
[building\_archetype](@ref) via the [building\_archetype\_\_building\_weather](@ref)
relationship. Regardless, the final [`WeatherData`](@ref) structs are stored
into a `weather_data_dictionary`, which can be examined through the Julia REPL
after the main program has finished.


## Process `ArchetypeBuilding` structs

With all the pieces now in place, the main program can finally process all the
data into lumped-capacitance thermal models depicting the desired synthetic
average archetype buildings. This is handled by the [`ArchetypeBuilding`](@ref)
struct constructor, and takes as input [The `building_archetype` definition](@ref),
as well as the appropriate [`ScopeData`](@ref) and [`WeatherData`](@ref)
processed during the previous steps.

The [`ArchetypeBuilding`](@ref) contains all the information about the final
lumped-capacitance thermal model of the synthetic average archetype building,
as well as the definitions used in its construction. The final
[`ArchetypeBuilding`](@ref) structs are stored into an `archetype_dictionary`,
which can be examined through the Julia REPL after the main program has finished.


## Solve the HVAC demand

After processing the [`ArchetypeBuilding`](@ref)s, the main program will
calculate a baseline/reference heating/cooling demand using the
lumped-capacitance thermal models of the synthetic average archetype buildings.
This is handled by the [`ArchetypeBuildingResults`](@ref) struct and its
constructor.

The final [`ArchetypeBuildingResults`](@ref) are written back into the
input *Datastore*, as well as stored into an `archetype_results_dictionary`,
which can be examined through the Julia REPL after the main program has finished.


## Export SpineOpt input data

If the `-spineopt <url>` argument is given, the main program will attempt to
convert the [`ArchetypeBuilding`](@ref)s in the `archetype_dictionary`
into *SpineOpt* energy system model input data, and export that input data into
the *Spine Datastore* at the given `url`.

The input data creation is handled by the [`SpineOptInput`](@ref) struct and
its constructor.


## Export Backbone input data

If the `-backbone <url>` argument is given, the main program will attempt to
convert the [`ArchetypeBuilding`](@ref)s in the `archetype_dictionary`
into *Backbone* energy system model input data, and export that input data into
the *Spine Datastore* at the given `url`.

The input data creation is handled by the [`BackboneInput`](@ref) struct and
its constructor.