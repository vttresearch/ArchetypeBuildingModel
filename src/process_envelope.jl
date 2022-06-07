#=
    process_envelope.jl

This file contains the functions for processing the building envelope.
=#


"""
    process_building_envelope(archetype::Object, data::ScopeData; mod::Module = Main)

Calculate the dimensions of the `archetype` building envelope based on the assumptions and
aggregated building stock data.

NOTE! The `mod` keyword changes from which Module data is accessed from,
`Main` by default.

Essentially, the archetype buildings are assumed to be rectangular in shape, with their external dimensions
primarily governed by the `number_of_storeys`, `building_frame_depth_m`, and `room_height_m` parameters.
The `average_gross_floor_area_m2_per_building` is divided into the desired number of floors,
and the perimeter walls are calculated by assuming `building_frame_depth_m` as the length of two parallel walls,
and determining the length of the two perpendicular walls based on the floor area.
In case of non-integer `number_of_storeys`, only the topmost floor is assumed to deviate in shape from the rest,
and the size of the partial floor is reduced by decreasing the length of the wall
that isn't fixed to `building_frame_depth_m`.

**Note that for interior structures like separating floors and partition walls,
the `surface_area_m2` only includes the surface area of one side of the structure!**
Accounting for coupling of both surfaces to the interior air is handled through the
`external_U_value_to_ambient_air_W_m2K` and `internal_U_value_to_structure_W_m2K` parameters.

The building envelope calculations proceed as follows:
1. Base floor dimensions using [`calculate_base_floor_dimensions`](@ref).
2. Roof dimensions using [`calculate_roof_dimensions`](@ref).
3. Separating floor dimensions using [`calculate_separating_floor_dimensions`](@ref).
4. Total vertical envelope area using [`calculate_vertical_envelope_surface_area`](@ref).
5. Window dimensions using [`calculate_window_dimensions`](@ref).
6. Exterior wall dimensions using [`calculate_exterior_wall_dimensions`](@ref).
7. Partition wall dimensions using [`calculate_partition_wall_dimensions`](@ref).
"""
function process_building_envelope(archetype::Object, data::ScopeData; mod::Module = Main)
    # Calculate envelope dimensions
    base_floor = calculate_base_floor_dimensions(
        data,
        mod.number_of_storeys(building_archetype = archetype),
        mod.building_frame_depth_m(building_archetype = archetype),
    )
    roof = calculate_roof_dimensions(
        base_floor,
        mod.number_of_storeys(building_archetype = archetype),
        mod.building_frame_depth_m(building_archetype = archetype),
    )
    separating_floor = calculate_separating_floor_dimensions(
        data,
        base_floor,
        mod.number_of_storeys(building_archetype = archetype),
        mod.building_frame_depth_m(building_archetype = archetype),
    )
    vertical_envelope_area_m2 = calculate_vertical_envelope_surface_area(
        base_floor,
        separating_floor,
        mod.room_height_m(building_archetype = archetype),
    )
    window = calculate_window_dimensions(
        vertical_envelope_area_m2,
        mod.window_area_to_external_wall_ratio_m2_m2(building_archetype = archetype),
    )
    exterior_wall, light_exterior_wall = calculate_exterior_wall_dimensions(
        vertical_envelope_area_m2,
        window,
        mod.external_wall_load_bearing_fraction(building_archetype = archetype),
        mod.number_of_storeys(building_archetype = archetype),
        mod.room_height_m(building_archetype = archetype),
    )
    partition_wall, light_partition_wall = calculate_partition_wall_dimensions(
        vertical_envelope_area_m2,
        mod.partition_wall_length_ratio_to_external_walls_m_m(
            building_archetype = archetype,
        ),
        mod.partition_wall_load_bearing_fraction(building_archetype = archetype),
    )

    # Return the envelope dimensions
    return base_floor,
    exterior_wall,
    light_exterior_wall,
    light_partition_wall,
    partition_wall,
    roof,
    separating_floor,
    window
end


"""
    calculate_base_floor_dimensions(
        data::ScopeData,
        storeys::Real,
        frame_depth_m::Real
    )

Calculate the base floor dimensions, assuming a rectangular building.

Essentially, calculates the surface area of the base floor in [m2],
as well as the thermal bridge length [m] based on the
known gross-floor area (GFA), the assumed number of storeys,
and the assumed frame depth of the building.
The thermal bridge length is simply assumed to correspond to the
length of the perimeter of the base floor.

```math
A_\\text{bf} = \\frac{A_\\text{GFA}}{n_\\text{storeys}}, \\\\
l_\\text{bf} = 2 \\left( \\frac{A_\\text{bf}}{d_\\text{frame}} + d_\\text{frame} \\right)
```
"""
function calculate_base_floor_dimensions(
    data::ScopeData,
    storeys::Real,
    frame_depth_m::Real,
)
    A_bf = data.average_gross_floor_area_m2_per_building / storeys
    return (
        linear_thermal_bridge_length_m = 2 * (A_bf / frame_depth_m + frame_depth_m),
        surface_area_m2 = A_bf,
    )
end


"""
    calculate_roof_dimensions(
        base_floor::NamedTuple,
        storeys::Real,
        frame_depth_m::Real
    )

Calculate the roof dimensions assuming a rectangular building.

Essentially, calculates the surface area of the roof in [m2],
as well as the thermal bridge length [m] based on the
dimensions of the base floor ([`calculate_base_floor_dimensions`](@ref)).
The thermal bridge length is simply assumed to correspond to the
length of the perimeter of the roof, and in case of a partial top floor,
the roof is divided into two separate surfaces,
thus increasing the thermal bridge length.

```math
A_\\text{r} = A_\\text{bf}, \\\\
l_\\text{r} = \\begin{cases}
l_\\text{bf}, \\qquad n_\\text{storey} \\in \\mathbb{N} \\\\
l_\\text{bf} + 2 d_\\text{frame}, \\qquad n_\\text{storey} \\notin \\mathbb{N}
\\end{cases}
```
"""
function calculate_roof_dimensions(
    base_floor::NamedTuple,
    frame_depth_m::Real,
    storeys::Real,
)
    return (
        linear_thermal_bridge_length_m = base_floor.linear_thermal_bridge_length_m +
                                         (mod(storeys, 1) != 0 && 2 * frame_depth_m),
        surface_area_m2 = base_floor.surface_area_m2,
    )
end


"""
    calculate_separating_floor_dimensions(
        data::ScopeData,
        base_floor::NamedTuple,
        storeys::Real,
        frame_depth_m::Real
    )

Calculate the separating floor dimensions assuming a rectangular building.

Essentially, calculates the *one-sided* surface area of the separating floors in [m2],
as well as the thermal bridge length [m] based on the known gross-floor area (GFA),
the assumed number of storeys, and the dimensions of the base floor ([`calculate_base_floor_dimensions`](@ref)).
The thermal bridge length is assumed to correspond to the
length of the perimeter of the separating floors.

```math
A_\\text{sf} = A_\\text{GFA} - A_\\text{bf}, \\\\
l_\\text{sf} = \\begin{cases}
\\left( \\lfloor n_\\text{storeys} \\rfloor - 1 \\right) l_\\text{bf}, \\quad & n_\\text{storeys} \\in \\mathbb{N} \\\\
\\left( \\lfloor n_\\text{storeys} \\rfloor - 1 \\right) l_\\text{bf} + 2 \\left( \\frac{A_\\text{bf} (n_\\text{storeys} - \\lfloor n_\\text{storeys} \\rfloor)}{d_\\text{frame}} + d_\\text{frame} \\right), \\quad & n_\\text{storeys} \\notin \\mathbb{N}
\\end{cases}
```
"""
function calculate_separating_floor_dimensions(
    data::ScopeData,
    base_floor::NamedTuple,
    storeys::Real,
    frame_depth_m::Real,
)
    full_storeys = floor(storeys)
    partial_storey = mod(storeys, 1)
    return (
        linear_thermal_bridge_length_m = (full_storeys - 1) *
                                         base_floor.linear_thermal_bridge_length_m +
                                         2 * (
            partial_storey != 0 &&
            (base_floor.surface_area_m2 * partial_storey / frame_depth_m + frame_depth_m)
        ),
        surface_area_m2 = data.average_gross_floor_area_m2_per_building -
                          base_floor.surface_area_m2,
    )
end


"""
    calculate_vertical_envelope_surface_area(
        base_floor::NamedTuple,
        separating_floor::NamedTuple,
        room_height_m::Real,
    )

Calculate the vertical envelope surface area [m2] assuming a rectangular building.

An auxiliary function used for window and exterior wall calculations.
```math
A_\\text{vertical envelope} = h_\\text{room} \\left( l_\\text{bf} + l_\\text{sf} \\right)
```
"""
function calculate_vertical_envelope_surface_area(
    base_floor::NamedTuple,
    separating_floor::NamedTuple,
    room_height_m::Real,
)
    return room_height_m * (
        base_floor.linear_thermal_bridge_length_m +
        separating_floor.linear_thermal_bridge_length_m
    )
end


"""
    calculate_window_dimensions(
        vertical_envelope_area_m2::Real,
        window_to_wall_ratio::Real
    )

Calculate the window dimensions assuming a rectangular building.

Calculates the window surface area [m2] simply using the assumed `window_to_wall_ratio`.
The linear thermal bridge lenght is set to zero, as it's not really applicable
without significantly better information about the number and size of the
individual windows.
```math
A_\\text{w} = w A_\\text{vertical envelope}, \\\\
l_\\text{w} = 0
```
"""
function calculate_window_dimensions(
    vertical_envelope_area_m2::Real,
    window_to_wall_ratio::Real,
)
    return (
        linear_thermal_bridge_length_m = 0.0,
        surface_area_m2 = vertical_envelope_area_m2 * window_to_wall_ratio,
    )
end


"""
    calculate_exterior_wall_dimensions(
        vertical_envelope_area_m2::Real,
        window::NamedTuple,
        load_bearing_fraction::Real,
        storeys::Real,
        room_height_m::Real
    )

Calculate the dimensions of exterior walls, assuming a rectangular building.

The exterior wall surface area [m2] is calculated simply as the difference between
`vertical_envelope_area_m2` ([`calculate_vertical_envelope_surface_area`](@ref))
and the window surface area ([`calculate_window_dimensions`](@ref)).
The linear thermal bridge length [m] includes the corners of the exterior walls.
The results are provided for load-bearing and light exterior walls separately,
based on the given `load_bearing_fraction`.
```math
A_\\text{ewlb} = f_\\text{lb} \\left( A_\\text{vertical envelope} - A_\\text{w} \\right), \\\\
l_\\text{ewlb} = 4 f_\\text{lb} \\lceil n_\\text{storeys} \\rceil h_\\text{room}
```
The dimensions of light exterior walls are calculated similarly to the above,
except that `(1-f_lb)` is used as the coefficient instead.

Returns the load-bearing exterior wall dimensions first,
and the light exterior wall dimensions second.
"""
function calculate_exterior_wall_dimensions(
    vertical_envelope_area_m2::Real,
    window::NamedTuple,
    load_bearing_fraction::Real,
    storeys::Real,
    room_height_m::Real,
)
    area_m2 = vertical_envelope_area_m2 - window.surface_area_m2
    thermal_bridge_length_m = 4 * ceil(storeys) * room_height_m
    return (
        linear_thermal_bridge_length_m = load_bearing_fraction * thermal_bridge_length_m,
        surface_area_m2 = load_bearing_fraction * area_m2,
    ),
    (
        linear_thermal_bridge_length_m = (1 - load_bearing_fraction) *
                                         thermal_bridge_length_m,
        surface_area_m2 = (1 - load_bearing_fraction) * area_m2,
    )
end


"""
    calculate_partition_wall_dimensions(
        vertical_envelope_area_m2::Real,
        partition_wall_ratio::Real,
        load_bearing_fraction::Real
    )

Calculate the partition wall dimensions, assuming a rectangular building.

The *one-sided* partition wall surface area [m2] is calculated based on the
`vertical_envelope_area_m2` ([`calculate_vertical_envelope_surface_area`](@ref)),
and the given `partition_wall_ratio`.
The linear thermal bridge length [m] is set to zero, as producing any actual
estimation would require a lot more information about the layout of the building.
The results are provided for load-bearing and light exterior walls separately,
based on the given `load_bearing_fraction`.
```math
A_\\text{pwlb} = f_\\text{lb} r_\\text{pw} A_\\text{vertical envelope}, \\\\
l_\\text{pwlb} = 0
```
The dimensions of light partition walls are calculated similarly to the above,
except that `(1-f_lb)` is used as the coefficient instead.

Returns the load-bearing partition wall dimensions first,
and the light partition wall dimensions second.
"""
function calculate_partition_wall_dimensions(
    vertical_envelope_area_m2::Real,
    partition_wall_ratio::Real,
    load_bearing_fraction::Real,
)
    area_m2 = partition_wall_ratio * vertical_envelope_area_m2
    return (
        linear_thermal_bridge_length_m = 0.0,
        surface_area_m2 = load_bearing_fraction * area_m2,
    ),
    (
        linear_thermal_bridge_length_m = 0.0,
        surface_area_m2 = (1 - load_bearing_fraction) * area_m2,
    )
end
