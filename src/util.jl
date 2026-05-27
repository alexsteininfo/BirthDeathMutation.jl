function get_simulation(input::SimulationInput, rng::AbstractRNG=Random.GLOBAL_RNG)
    return runsimulation(input, rng)
end
