struct SampledData
    VAF::Vector{Float64}
    depth::Vector{Int64}
end

abstract type SimulationResult end

struct Simulation{S<:SimulationInput, T<:AbstractPopulation} <: SimulationResult
    input::S
    output::T
end

get_simulation(sim::Simulation) = sim

function Base.show(io::IO, simulation::SimulationResult)
    @printf(io, "===================================================================\n")
    show(io, simulation.input)
    @printf(io, "===================================================================\n")
    show(io, simulation.output)
    @printf(io, "\n===================================================================")
end
