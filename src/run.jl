"""
    runsimulation(input, [selection], rng; timefunc=exptime, returnextinct=false)

Run a single simulation defined by `input` and `selection` (defaults to
`NeutralSelection()`).
"""
function runsimulation end

function runsimulation(input::SimulationInput, rng::AbstractRNG=Random.GLOBAL_RNG; kwargs...)
    return runsimulation(input, NeutralSelection(), rng; kwargs...)
end

function runsimulation(
    input::SimulationInput,
    selection::AbstractSelection,
    rng::AbstractRNG=Random.GLOBAL_RNG;
    timefunc=exptime,
    returnextinct=false
)
    population = initialize_population(input)
    while true
        nextID = initialize_counters(population)
        population, = simulate!(population, input, selection, nextID, rng; timefunc)
        if length(allcells(population)) != 0 || returnextinct
            break
        else
            population = initialize_population(input)
        end
    end
    return Simulation(input, population)
end

"""
    runsimulation_timeseries_returnfinalpop(input, [selection], timesteps, func, rng;
        timefunc=exptime, returnextinct=false)

Run a single simulation. Calls `func(population)` at each timestep and returns results
as a vector along with the final population state.
"""
function runsimulation_timeseries_returnfinalpop end

function runsimulation_timeseries_returnfinalpop(
    input::SimulationInput,
    timesteps,
    func,
    rng::AbstractRNG=Random.GLOBAL_RNG;
    kwargs...
)
    return runsimulation_timeseries_returnfinalpop(
        input, NeutralSelection(), timesteps, func, rng; kwargs...
    )
end

function runsimulation_timeseries_returnfinalpop(
    input::SimulationInput,
    selection::AbstractSelection,
    timesteps,
    func,
    rng::AbstractRNG=Random.GLOBAL_RNG;
    timefunc=exptime,
    returnextinct=false
)
    population = initialize_population(input)
    data = []
    while true
        nextID = initialize_counters(population)
        t0 = 0.0
        for t in timesteps
            simulate!(population, input, selection, nextID, rng; timefunc, t0, tmax=t)
            push!(data, func(population))
            t0 = t
        end
        if length(allcells(population)) != 0 || returnextinct
            break
        else
            population = initialize_population(input)
            data = []
        end
    end
    return data, Simulation(input, population)
end

"""
    runsimulation_timeseries(input, [selection], timesteps, func, rng;
        timefunc=exptime, returnextinct=false)

Run a single simulation. Calls `func(population)` at each timestep and returns results
as a vector.
"""
function runsimulation_timeseries end

function runsimulation_timeseries(
    input::SimulationInput,
    timesteps,
    func,
    rng::AbstractRNG=Random.GLOBAL_RNG;
    kwargs...
)
    return runsimulation_timeseries(
        input, NeutralSelection(), timesteps, func, rng; kwargs...
    )
end

function runsimulation_timeseries(
    input::SimulationInput,
    selection::AbstractSelection,
    timesteps,
    func,
    rng::AbstractRNG=Random.GLOBAL_RNG;
    kwargs...
)
    return runsimulation_timeseries_returnfinalpop(
        input, selection, timesteps, func, rng; kwargs...
    )[1]
end

getinputrates(input::BranchingInput) = input.birthrate, input.deathrate, 0.0, 0.0
getinputrates(input::MoranInput) = 0.0, 0.0, input.moranrate, 0.0
getinputrates(input::BranchingMoranInput) = input.birthrate, input.deathrate, input.moranrate, 0.0
