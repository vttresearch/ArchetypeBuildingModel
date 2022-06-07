#=
    solve_demand.jl

Functions for solving the heating/cooling demands of `ArchetypeBuilding`s.
=#


"""
    solve_heating_demand(archetype::ArchetypeBuilding, free_dynamics::Bool)

Solve the heating/cooling demand of the `archetype`.

Note that this function calculates the "final energy demand" of the archetype
building, and not the energy consumption of it's HVAC systems.
See the [`solve_consumption`](@ref) function for that.
Essentially, performs the following steps:
1. Check external load data and [`determine_temporal_structure`](@ref).
2. [`form_and_invert_dynamics_matrix`](@ref) for the free temperature dynamics.
3. Initialize external load an thermal mass vectors using [`initialize_rhs`](@ref).
4. Initialize temperature and temperature limit vectors using [`initialize_temperatures`](@ref).
5. Solve the heating demand using [`solve_heating_demand_loop`](@ref).
6. Rearrange the solution into `Dict`s and return them.

Uses an extremely simple rule-based control to solve the heating/cooling
demand of the archetype building. The controller intervenes
whenever node temperatures would deviate from permitted limits,
and provides the required energy input to maintain the system at the limit.

The building dynamics are discretized using implicit *(backwards)* Euler,
mainly for consistency with our existing energy system modelling tools
like Backbone or SpineOpt. In principle, I believe the system could be solved
analytically similar to my Master's Thesis:
*Energy Management in Households with Coupled Photovoltaics and Electric Vehicles, Topi Rasku, 2015, Aalto University School of Science*.
"""
function solve_heating_demand(archetype::ArchetypeBuilding, free_dynamics::Bool)
    # Check that the external load data for the abstract nodes makes sense,
    # and determine the temporal scope and resolution for the simulation.
    indices, delta_t = determine_temporal_structure(archetype)

    # Form and invert the implicit Euler free dynamics matrix.
    dynamics_matrix, inverted_dynamics_matrix =
        form_and_invert_dynamics_matrix(archetype, delta_t)

    # Initialize the `external_load` and thermal mass vectors.
    external_load_vector, thermal_mass_vector = initialize_rhs(archetype, indices, delta_t)

    # Initialize the temperature vector and the temperature limit vectors.
    initial_temperatures, min_temperatures, max_temperatures = initialize_temperatures(
        archetype,
        indices,
        dynamics_matrix,
        inverted_dynamics_matrix,
        external_load_vector,
        thermal_mass_vector,
        free_dynamics,
    )

    # Solve the heating demand for the entire set of indices.
    temperatures, hvac_demand = solve_heating_demand_loop(
        indices,
        dynamics_matrix,
        inverted_dynamics_matrix,
        initial_temperatures,
        min_temperatures,
        max_temperatures,
        external_load_vector,
        thermal_mass_vector,
        free_dynamics,
    )

    # Rearrange the results into Dicts,
    # and remove initial values from the temperature array.
    nodes = keys(archetype.abstract_nodes)
    init_temp_dict = Dict(zip(nodes, popfirst!(temperatures)))
    temp_dict = Dict(
        zip(
            nodes,
            [
                TimeSeries(indices, getindex.(temperatures, i), false, false) for
                (i, n) in enumerate(nodes)
            ],
        ),
    )
    hvac_demand_dict = Dict(
        zip(
            nodes,
            [
                TimeSeries(indices, getindex.(hvac_demand, i), false, false) for
                (i, n) in enumerate(nodes)
            ],
        ),
    )

    # Return the solution.
    return init_temp_dict, temp_dict, hvac_demand_dict
end


"""
    determine_temporal_structure(archetype::ArchetypeBuilding)

Check that `external_load` timeseries are consistent in the `AbstractNodeNetwork`,
and determine the time series indices and the `delta_t`.

Note that the time series need to have a constant `delta_t` in order for the
dynamic matrix to be time-invarying, speeding up the solving process significantly.
"""
function determine_temporal_structure(archetype::ArchetypeBuilding)
    # Check that all nodes have identical `external_load` time series indices.
    indices = first(archetype.abstract_nodes)[2].external_load.indexes
    if !all(n.external_load.indexes == indices for (k, n) in archetype.abstract_nodes)
        return @error """
        `external_load` time series are indexed different for `abstract_nodes`
        of `archetype_building` `$(archetype)`!
        """
    end

    # Calculate the delta t in hours, all time steps need to have constant value.
    delta_t = getfield.(Hour.(diff(indices)), :value)
    if !all(delta_t .== first(delta_t))
        return @error """
        `external_load` time series of `archetype` `$(archetype)` must have
        a constant time step length!
        """
    end

    return indices, first(delta_t)
end


"""
    form_and_invert_dynamics_matrix(archetype::ArchetypeBuilding, delta_t::Int64)

Forms and inverts the implicit Euler discretized dynamics matrix for the
`AbstractNodeNetwork`.
"""
function form_and_invert_dynamics_matrix(archetype::ArchetypeBuilding, delta_t::Int64)
    # Initialize the implicit Euler dynamics matrix.
    len = length(archetype.abstract_nodes)
    M = Matrix{Float64}(undef, len, len)

    # Loop over the matrix to fill in the values.
    enum_nodes = enumerate(archetype.abstract_nodes)
    for (i, (k1, n1)) in enum_nodes
        for (j, (k2, n2)) in enum_nodes
            if i == j
                M[i, j] =
                    n1.thermal_mass_Wh_K / delta_t +
                    n1.self_discharge_coefficient_W_K +
                    reduce(+, values(n1.heat_transfer_coefficients_W_K); init = 0.0)
            else
                M[i, j] = -get(n1.heat_transfer_coefficients_W_K, k2, 0.0)
            end
        end
    end

    # Return the dynamics matrix and its inverse
    return M, inv(M)
end


"""
    initialize_temperatures(
        archetype::ArchetypeBuilding,
        indices::Vector{DateTime},
        dynamics_matrix::Matrix{Float64},
        inverted_dynamics_matrix::Matrix{Float64},
        external_load_vector::Vector{Vector{Float64}},
        thermal_mass_vector::Vector{Float64},
        free_dynamics::Bool,
    )

Initialize the temperature and temperature limit vectors for the heating/cooling
demand calculations.

Initial temperatures are solved by repeatedly solving the first 24 hours
until the end-result no longer changes.
The initialization is abandoned if no stable initial temperatures are found
within a thousand 24-hour solves.
In this case, the minimum permitted temperatures are used as the initial
temperatures for each node.
"""
function initialize_temperatures(
    archetype::ArchetypeBuilding,
    indices::Vector{DateTime},
    dynamics_matrix::Matrix{Float64},
    inverted_dynamics_matrix::Matrix{Float64},
    external_load_vector::Vector{Vector{Float64}},
    thermal_mass_vector::Vector{Float64},
    free_dynamics::Bool,
)
    # Initialize the temperature vector using minimum allowed temperatures.
    initial_temperatures =
        float.([n.minimum_temperature_K for (k, n) in archetype.abstract_nodes])

    # Form temperature limit vectors for convenience.
    min_temperatures = deepcopy(initial_temperatures)
    max_temperatures =
        float.([n.maximum_temperature_K for (k, n) in archetype.abstract_nodes])

    # Solve the initial temperatures via repeatedly solving the first 24 hours
    # until the temperatures converge, starting from the permitted minimums.
    for i = 1:1000
        temps, hvac = solve_heating_demand_loop(
            indices[1:24],
            dynamics_matrix,
            inverted_dynamics_matrix,
            initial_temperatures,
            min_temperatures,
            max_temperatures,
            external_load_vector,
            thermal_mass_vector,
            free_dynamics,
        )
        if isapprox(last(temps), initial_temperatures)
            println("Stable initial temperatures found with $(i) iterations.")
            return last(temps), min_temperatures, max_temperatures
        end
        initial_temperatures = last(temps)
    end
    println("""
            No stable initial temperatures found!
            Using minimum permitted temperatures instead.
            """)
    return deepcopy(min_temperatures), min_temperatures, max_temperatures
end


"""
    initialize_rhs(
        archetype::ArchetypeBuilding,
        indices::Vector{Dates.DateTime},
        delta_t::Int64
    )

Initialize the right-hand side of the linear equation system,
meaning the impact of the `external_load` and previous temperatures.
"""
function initialize_rhs(
    archetype::ArchetypeBuilding,
    indices::Vector{Dates.DateTime},
    delta_t::Int64,
)
    # Process the nodal `external_loads` into a nested vector for easy access.
    external_load_vector = [
        [n.external_load.values[i] for (k, n) in archetype.abstract_nodes] for
        (i, t) in enumerate(indices)
    ]

    # Calculate the thermal mass vector to account for previous temperatures.
    thermal_mass_vector =
        [n.thermal_mass_Wh_K / delta_t for (k, n) in archetype.abstract_nodes]

    return external_load_vector, thermal_mass_vector
end


"""
    solve_heating_demand_loop(
        indices::Vector{DateTime},
        dynamics_matrix::Matrix{Float64},
        inverted_dynamics_matrix::Matrix{Float64},
        temperatures::Vector{Vector{Float64}},
        min_temperatures::Vector{Float64},
        max_temperatures::Vector{Float64},
        external_load_vector::Vector{Vector{Float64}},
        thermal_mass_vector::Vector{Float64},
        free_dynamics::Bool,
        initial_temperatures::Vector{Float64}
    )

Solve the heating/cooling demand one timestep at a time over the given indices.

Essentially, performs the following steps:
1. Initialize the temperature vector, HVAC demand vector, and a dictionary for the dynamic matrices for solving the problem.
2. Loop over the given `indices` and do the following:
    3. Solve new temperatures if HVAC not in use.
    4. Check if new temperatures would violate temperature limits.
    5. If necessary, solve the HVAC demand required to keep temperatures within set limits.
6. Return the solved temperatures and HVAC demand for each node and index.
"""
function solve_heating_demand_loop(
    indices::Vector{DateTime},
    dynamics_matrix::Matrix{Float64},
    inverted_dynamics_matrix::Matrix{Float64},
    initial_temperatures::Vector{Float64},
    min_temperatures::Vector{Float64},
    max_temperatures::Vector{Float64},
    external_load_vector::Vector{Vector{Float64}},
    thermal_mass_vector::Vector{Float64},
    free_dynamics::Bool,
)
    # Initialize the temperature vector, the HVAC heating/cooling demand vector,
    # and the heating/cooling demand solving matrix Dict
    temperatures = [initial_temperatures]
    hvac_demand = Vector{Vector{Float64}}()
    inverse_hvac_matrices = Dict{Vector{Int64},Matrix{Float64}}()

    # Loop over the indices, and solve the dynamics/HVAC demand.
    for (i, t) in enumerate(indices)
        # Calculate the new temperatures without HVAC.
        previous_temperature_effect_vector = thermal_mass_vector .* last(temperatures)
        new_temperatures =
            inverted_dynamics_matrix *
            (external_load_vector[i] + previous_temperature_effect_vector)

        # Check if the temperatures are within permissible limits.
        max_temp_check = new_temperatures .<= max_temperatures
        min_temp_check = new_temperatures .>= min_temperatures
        temp_check = max_temp_check .* min_temp_check
        if free_dynamics || all(temp_check)
            # If yes, simply save the new temperatures & zero demand, and move on.
            push!(temperatures, new_temperatures)
            push!(hvac_demand, zeros(size(new_temperatures)))
        else
            # Else, calculate the required HVAC demand.
            # Check if we've already calculated the inverse HVAC matrix for this case,
            # if not, calculate and store it for future reference.
            inverse_hvac_matrix = get(inverse_hvac_matrices, temp_check, nothing)
            if isnothing(inverse_hvac_matrix)
                inverse_hvac_matrix =
                    form_and_invert_hvac_matrix(dynamics_matrix, temp_check)
                inverse_hvac_matrices[temp_check] = inverse_hvac_matrix
            end

            # Find which nodes violate the temperature limits.
            fixed_max_temp_inds = findall(.!(max_temp_check))
            fixed_min_temp_inds = findall(.!(min_temp_check))

            # Solve the HVAC demand and new temperatures.
            # This is a bit complicated due to the violated temperatures
            # becoming fixed, moving them from the left to the right-hand side
            # of the system of equations.
            hvac_solution =
                inverse_hvac_matrix * (
                    external_load_vector[i] .+ previous_temperature_effect_vector .-
                    reduce(
                        .+,
                        dynamics_matrix[:, j] .* min_temperatures[j] for
                        j in fixed_min_temp_inds;
                        init = 0.0,
                    ) .- reduce(
                        .+,
                        dynamics_matrix[:, j] .* max_temperatures[j] for
                        j in fixed_max_temp_inds;
                        init = 0.0,
                    )
                )

            # Finally, separate the results into temperature and HVAC vectors.
            # Note that violated temperatures were fixed, and replaced with
            # the HVAC demand variables.
            new_temperatures = deepcopy(hvac_solution)
            new_temperatures[fixed_max_temp_inds] = max_temperatures[fixed_max_temp_inds]
            new_temperatures[fixed_min_temp_inds] = min_temperatures[fixed_min_temp_inds]
            hvac = zeros(size(hvac_solution))
            hvac[fixed_max_temp_inds] = hvac_solution[fixed_max_temp_inds]
            hvac[fixed_min_temp_inds] = hvac_solution[fixed_min_temp_inds]
            push!(temperatures, new_temperatures)
            push!(hvac_demand, hvac)
        end
    end
    return temperatures, hvac_demand
end


"""
    form_and_invert_hvac_matrix(
        dynamics_matrix::Matrix{Float64},
        temp_check::BitVector
    )

Forms and inverts the matrix for solving HVAC demand in different situations.
"""
function form_and_invert_hvac_matrix(
    dynamics_matrix::Matrix{Float64},
    temp_check::BitVector,
)
    # HVAC matrix is based on the dynamics matrix,
    # except that all violated temperatures are fixed and their variables
    # are replaced with hvac demand.
    hvac_matrix = deepcopy(dynamics_matrix)
    for i in findall(.!(temp_check))
        hvac_matrix[:, i] .= 0.0
        hvac_matrix[i, i] = -1.0
    end

    return inv(hvac_matrix)
end
