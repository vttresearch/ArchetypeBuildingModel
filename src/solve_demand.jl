#=
    solve_demand.jl

Functions for solving the heating/cooling demands of `ArchetypeBuilding`s.
=#


"""
    solve_heating_demand(
        archetype::ArchetypeBuilding,
        free_dynamics::Bool,
        initial_temperatures::Union{Nothing,Dict{Object,Float64}};
        realization::Symbol = :realization,
    )

Solve the heating/cooling demand of the `archetype`.

Note that this function calculates the "final energy demand" of the archetype
building, and not the energy consumption of it's HVAC systems.
Furthermore, the calculations are deterministic, with `realization` defining
the true data from potentially stochastic input.
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
analytically similar to my Master's Thesis, if better accuracy would be desired:
*Energy Management in Households with Coupled Photovoltaics and Electric Vehicles, Topi Rasku, 2015, Aalto University School of Science*.

The idea of solving heating demand calculations can be as follows,
starting from the energy balance equation for node `n`
```math
C_n \\frac{dT_n(t)}{dt} = - \\rho_n T_n(t) + \\sum_{m \\in N} \\left[ H_{n,m} \\left( T_m(t) - T_n(t) \\right) \\right] + \\Phi_n(t),
```
where `C_n` is the heat capacity of node `n`,
`T_n(t)` is the time-dependent temperature of node `n` at time `t`,
`ρ_n` is the self-discharge from node `n`,
`N` is the set of temperature nodes in the lumped-capacitance model of the building,
`H_n,m` is the heat transfer coefficient between nodes `n` and `m`,
and `Φ_n(t)` is the time-dependent total external heat load on node `n`.
Using the implicit Euler discretization, the above can be cast into
```math
\\left( \\frac{C_n}{\\Delta t} + \\rho_n + \\sum_{m \\in N} H_{n,m} \\right) T_{n,t} - \\sum_{m \\in N} \\left[ H_{n,m} T_{m,t} \\right] = \\Phi_{n,t} + \\frac{C_n}{\\Delta t} T_{n,t-\\Delta t},
```
where `Δt` is the length of the discretized time step.
Since we always know the temperatures on the previous time step `t-Δt`,
the above can be expressed in matrix form and solved as
```math
\\bm{A} \\hat{T} = \\hat{\\Phi}, \\\\
\\hat{T} = \\bm{A}^{-1} \\hat{\\Phi},
```
where `A` is the so-called *dynamics matrix*, 
`T` is the current temperature vector,
and `Φ` is the right-hand side vector, containing the effect of external loads
and previous temperatures.

The above is used to calculate the temperatures of the nodes on each subsequent
time step. However, when any of the node temperatures would violate the defined
minimum and maximum permitted temperatures, that temperature variable is instead
fixed to the violated boundary, moved to the right-hand side,
and replaced with a heating/cooling demand variable `ϕ_m` instead.
This results in a slightly modified problem to be solved
```math
\\hat{\\phi} = \\left( \\bm{A} - \\sum_{m \\in M}[\\bm{A}_{m} + \\bm{I}_{m}] \\right)^{-1} \\left( \\hat{\\Phi} - \\sum_{m \\in M}[\\hat{A}_{m} T'_m] \\right),
```
where `ϕ` is the modified temperature vector with the fixed temperature replaced
with the heating/cooling demand variable,
`m ∈ M` are the nodes that would violate their permitted bounds,
and are therefore fixed at the boundary `T'_m`,
`A_m` represents a selection of column `m` from the matrix `A`,
and `I_m` represents column `m` from an identity matrix.
Please note that there are both matrix and vector selections of `A_m`,
where the matrix selection preserves the dimensions with filling zeroes,
while the vector selection is essentially only the selected column in vector form.
"""
function solve_heating_demand(
    archetype::ArchetypeBuilding,
    free_dynamics::Bool,
    initial_temperatures::Union{Nothing,Dict{Object,Float64}};
    realization::Symbol=:realization
)
    # Check that the external load data for the abstract nodes makes sense,
    # and determine the temporal scope and resolution for the simulation.
    indices, delta_t = determine_temporal_structure(archetype; realization=realization)

    # Form and invert the implicit Euler free dynamics matrix.
    dynamics_matrix, inverted_dynamics_matrix =
        form_and_invert_dynamics_matrix(archetype, delta_t)

    # Initialize the `external_load` and thermal mass vectors.
    external_load_vector, thermal_mass_vector =
        initialize_rhs(archetype, indices, delta_t; realization=realization)

    # Initialize the temperature vector and the temperature limit vectors.
    init_temperatures, min_temperatures, max_temperatures = initialize_temperatures(
        archetype,
        indices,
        dynamics_matrix,
        inverted_dynamics_matrix,
        external_load_vector,
        thermal_mass_vector,
        free_dynamics,
        initial_temperatures,
    )

    # Solve the heating demand for the entire set of indices.
    temperatures, hvac_demand = solve_heating_demand_loop(
        indices,
        dynamics_matrix,
        inverted_dynamics_matrix,
        init_temperatures,
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
    determine_temporal_structure(
        archetype::ArchetypeBuilding;
        realization::Symbol = :realization,
    )

Check that `external_load` timeseries are consistent in the `AbstractNodeNetwork`,
and determine the time series indices and the `delta_t`.

Note that the time series need to have a constant `delta_t` in order for the
dynamic matrix to be time-invarying, speeding up the solving process significantly.
The `realization` keyword is necessary to indicate the true data from potentially
stochastic input.
"""
function determine_temporal_structure(
    archetype::ArchetypeBuilding;
    realization::Symbol=:realization
)
    # Check that all nodes have identical `external_load` time series indices.
    indices =
        parameter_value(first(archetype.abstract_nodes)[2].external_load)(
            scenario=realization,
        ).indexes
    if !all(
        parameter_value(n.external_load)(scenario=realization).indexes == indices for
        (k, n) in archetype.abstract_nodes
    )
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

The implicit Euler discretized dynamics matrix `A` is formed as follows:
```math
\\bm{A}_{n,m} = \\begin{cases}
\\frac{C_m}{\\Delta t} + \\rho_m + \\sum_{n' \\in N} H_{n,m}, \\qquad n = m, \\\\
- H_{n,m}, \\qquad n \\neq m,
\\end{cases}, \\quad \\text{where } n, m \\in N
```
where `A_n,m` is the element of the dynamic matrix `A` on row `n` and column `m`,
`C_m` is the thermal mass of node `m`,
`Δt` is the length of the discretized time step,
`ρ_m` is the self-discharge coefficient of node `m`,
`N` is the set of nodes included in the lumped-capacitance thermal model,
and `H_{m,n}` is the heat transfer coefficient between nodes `n` and `m`.
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
                    reduce(+, values(n1.heat_transfer_coefficients_W_K); init=0.0)
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
        initial_temperatures::Union{Nothing,Dict{Object,Float64}},
    )

Initialize the temperature and temperature limit vectors for the heating/cooling
demand calculations.

Initial temperatures are solved by repeatedly solving the first 24 hours
until the end-result no longer changes.
The initialization is abandoned if no stable initial temperatures are found
within a thousand 24-hour solves.
In this case, the minimum permitted temperatures are used as the initial
temperatures for each node, unless otherwise specified via `initial_temperatures`.
Internally, uses the [`solve_heating_demand_loop`](@ref) function.

See the [`solve_heating_demand`](@ref) function for the overall logic and
formulation of the heating demand calculations.
"""
function initialize_temperatures(
    archetype::ArchetypeBuilding,
    indices::Vector{DateTime},
    dynamics_matrix::Matrix{Float64},
    inverted_dynamics_matrix::Matrix{Float64},
    external_load_vector::Vector{Vector{Float64}},
    thermal_mass_vector::Vector{Float64},
    free_dynamics::Bool,
    initial_temperatures::Union{Nothing,Dict{Object,Float64}},
)
    # Fetch the allowed temperature limits.
    min_temperatures =
        float.([n.minimum_temperature_K for (k, n) in archetype.abstract_nodes])
    max_temperatures =
        float.([n.maximum_temperature_K for (k, n) in archetype.abstract_nodes])

    # Form the initial temperature vector.
    # Based on minimum temperatures unless otherwise defined.
    if !isnothing(initial_temperatures)
        init_temperatures =
            float.([
                get(initial_temperatures, n, min_temperatures[i]) for
                (i, n) in enumerate(keys(archetype.abstract_nodes))
            ])
        min_init_temperatures = deepcopy(init_temperatures)
        max_init_temperatures = deepcopy(max_temperatures)
        fixed_inds = findall(init_temperatures .!= min_temperatures)
        max_init_temperatures[fixed_inds] = init_temperatures[fixed_inds]
    else
        init_temperatures = deepcopy(min_temperatures)
        min_init_temperatures = deepcopy(min_temperatures)
        max_init_temperatures = deepcopy(max_temperatures)
    end

    # Solve the initial temperatures via repeatedly solving the first 24 hours
    # until the temperatures converge, starting from the permitted minimums.
    for i = 1:1000
        temps, hvac = solve_heating_demand_loop(
            indices[1:24],
            dynamics_matrix,
            inverted_dynamics_matrix,
            init_temperatures,
            min_init_temperatures,
            max_init_temperatures,
            external_load_vector,
            thermal_mass_vector,
            free_dynamics,
        )
        if isapprox(last(temps), init_temperatures)
            println("Stable initial temperatures found with $(i) iterations.")
            return last(temps), min_temperatures, max_temperatures
        end
        init_temperatures = last(temps)
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
        delta_t::Int64;
        realization::Symbol = :realization,
    )

Initialize the right-hand side of the linear equation system,
meaning the impact of the `external_load` and previous temperatures.

The `realization` keyword is used to indicate the true data from potentially
stochastic input.

See the [`solve_heating_demand`](@ref) function for the overall formulation.
This function returns the right-hand side components separately
```math
\\hat{\\Phi} = \\hat{\\Phi'} + \\hat{\\frac{C}{\\Delta t} T_{t-\\Delta t}},
```
where `Φ'` is the component of external loads,
and the rest is the component of the impact of previous temperatures.
The components are useful for the [`solve_heating_demand_loop`](@ref) function.
"""
function initialize_rhs(
    archetype::ArchetypeBuilding,
    indices::Vector{Dates.DateTime},
    delta_t::Int64;
    realization::Symbol=:realization
)
    # Process the nodal `external_loads` into a nested vector for easy access.
    external_load_vector = [
        [
            parameter_value(n.external_load)(scenario=realization).values[i] for
            (k, n) in archetype.abstract_nodes
        ] for (i, t) in enumerate(indices)
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

See the [`solve_heating_demand`](@ref) function for the overall formulation.
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
                        init=0.0
                    ) .- reduce(
                        .+,
                        dynamics_matrix[:, j] .* max_temperatures[j] for
                        j in fixed_max_temp_inds;
                        init=0.0
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

Essentially, this function performs the
```math
\\left( \\bm{A} - \\sum_{m \\in M}[\\bm{A}_{m} + \\bm{I}_{m}] \\right)^{-1}
```
transformation of the dynamics matrix `A`,
where the otherwise violated temperature variables `m ∈ M` are fixed
and replaced with a variable for the required heating/cooling demand.
See the [`solve_heating_demand`](@ref) function for the overall formulation.
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
