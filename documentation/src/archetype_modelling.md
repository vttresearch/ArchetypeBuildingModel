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


## Forming the building envelope


## Preparing the building loads


## Calculating the properties of the lumped-capacitance thermal nodes


## Calculating the properties of the HVAC equipment

