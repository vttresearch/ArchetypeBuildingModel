# Archetype building modelling

This section aims to explain the lumped-capacitance *(resistance-capacitance)*
modelling approach and the inherent assumptions used by this module in more detail,
aimed at a more expert audience interested in the inner workings of the module.
However, a lot of the techical details are not reproduced here,
but instead reference the relevant docstrings in the [Library](@ref) section.
Thus, this section still 

Overall, the structure of this section more or less follows the
[`ArchetypeBuilding`](@ref) constructor,
as it is the key `struct` containing and calculating all the information 
necessary for modelling an archetype building.
The [Process `ScopeData` structs](@ref), [Process `WeatherData` structs](@ref),
and [ArchetypeBuildingWeather.py](@ref) sections already sufficiently explain
the handling of the building stock statistics and weather data as input,
so we'll start with [Forming the building envelope](@ref) instead.
Next, we proceed to [Preparing the building loads](@ref),
namely the internal heat gains and solar heat gains for the archetype building.
Then, now that we know the dimensions of the building envelope and the
external loads on the building, we can proceed with
[Calculating the properties of the lumped-capacitance thermal nodes](@ref),
as well as with [Calculating the properties of the HVAC equipment](@ref).
Finally, `ArchetypeBuildingModel.jl` also includes a very simple rule-based
method for [Solving the heating demand and HVAC equipment consumption](@ref).
While the main goal of this module is to provide input data for optimization
models like [Backbone](https://cris.vtt.fi/en/publications/backbone) or
[SpineOpt](https://github.com/Spine-project/SpineOpt.jl),
having access to simple standalone baseline solutions is a quite helpful.


## Forming the building envelope

The first step in creating an archetype building is forming the building
envelope by calculating the dimensions of the different structures,
and storing the results in a corresponding [`EnvelopeData`](@ref) `struct`.
However, as `ArchetypeBuildingModel.jl` aims to remain useable on the building
stock scale, the geometry of the archetype buildings is heavily simplified.
The key assumptions for forming the archetype building envelope are as follows:

- The archetype buildings are assumed to be rectangular in shape, with an `average_gross_floor_area_m2_per_building` based on the appropriate [`ScopeData`](@ref).
- The key [building\_archetype](@ref) parameters controlling the shape of the building are:
    - [number\_of\_storeys](@ref): The assumed average number of storeys of the archetype building. In case of non-integer values, only the topmost separating floor is assumed to differ from the rest. Essentially, the `average_gross_floor_area_m2_per_building` is divided into this many floors.
    - [building\_frame\_depth\_m](@ref): The assumed average depth of the archetype building. Since dwellings typically have requirements for natural light *(at least in Finland)*, there's a limit for how deep the building frame can reasonably be when fenestrated from both sides. Naturally, this varies from building to building depending on their age and exact shape, but in general, I find it more reliable than using e.g. some assumption about the width-to-depth ratio of the buildings.
    - [room\_height\_m](@ref): The assumed average height of the rooms. In combination with the [number\_of\_storeys](@ref), this parameter determines how tall the archetype buildings are. For simplicity, all storeys are assumed to be equally high, and non-integer [number\_of\_storeys](@ref) are rounded up for estimating the height of the building.

See the [`ArchetypeBuildingModel.process_building_envelope`](@ref) docstring
for a detailed explanation of how the envelope shape related parameters
affect the calculations.
The actual equations for the surface areas and linear thermal
bridge lengths of the structures can be found from the docstrings of the
dedicated functions linked therein.


## Preparing the building loads


## Calculating the properties of the lumped-capacitance thermal nodes


## Calculating the properties of the HVAC equipment


## Solving the heating demand and HVAC equipment consumption