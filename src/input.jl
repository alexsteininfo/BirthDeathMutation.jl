"""
    SimulationInput

Supertype for simulation inputs.

# Subtypes
- `BranchingInput <: SinglelevelInput`: input for simulating a branching process
- `MoranInput <: SinglelevelInput`: input for simulating a Moran process
- `BranchingMoranInput <: SinglelevelInput`: input for simulating a branching process
    until fixed size is reached, then switching to Moran process
"""
abstract type SimulationInput end
abstract type SinglelevelInput <: SimulationInput end

"""
    BranchingInput <: SinglelevelInput <:SimulationInput

Input type for a single level branching process simulation that starts with a single cell.

# Fields:
- `Nmax::Int64 = 1000`: maximum number of cells
- `tmax::Float64 = Inf`: maximum time to run simulation
- `birthrate::Float64 = 1.0`: birth rate for wild-type cells
- `deathrate::Float64 = 0.0`: death rate for wild-type cells
- `clonalmutations::Int64 = 0`: number of mutations shared by all cells
- `μ::Vector{Float64} = [1.0]`: mutation rate per division per cell. Can be passed as a
    single `Float64`, if there is only one mutational process. Multiple values indicate
    multiple simulataneous processes.
- `mutationdist::Vector{Symbol} = [:poisson]`: defines the distibution for new
    mutations (:poisson, :fixed, :poissontimedep, :fixedtimedep, :geometric). Length should
    match `length(μ)`.
- `ploidy::Int64 = 2`: cell ploidy
"""
struct BranchingInput <: SinglelevelInput
    Nmax::Int64
    tmax::Float64
    birthrate::Float64
    deathrate::Float64
    clonalmutations::Int64
    μ::Vector{Float64}
    mutationdist::Vector{Symbol}
    ploidy::Int64
end

function BranchingInput(;
    Nmax = 1000,
    tmax = Inf,
    birthrate = 1.0,
    deathrate = 0.0,
    clonalmutations = 0,
    μ = [1.0],
    mutationdist = fill(:poisson, length(μ)),
    ploidy = 2
)
    μ = tovector(μ)
    mutationdist = tovector(mutationdist)
    @assert length(μ) == length(mutationdist) "μ and mutationdist are not same length"
    return BranchingInput(
        Nmax, tmax, birthrate, deathrate, clonalmutations, μ, mutationdist, ploidy
    )
end

tovector(a::T) where T = T[a]
tovector(a::Vector{T}) where T = a

"""
    MoranInput <: SinglelevelInput

Input type for a single level Moran process simulation that starts with `N` identical cells.

# Fields:
- `N::Int64 = 1000`: number of cells
- `tmax::Float64 = 15.0`: maximum time to run simulation
- `moranrate::Float64 = 1.0`: moran update rate for wild-type cells
- `moranincludeself::Bool = true`: determines whether the same cell can be chosen to both
    divide and die in a moran step (in which case one offspring is killed)
- `clonalmutations::Int64 = 0`: number of mutations shared by all cells
- `μ::Vector{Float64} = [1.0]`: mutation rate per division per cell. Can be passed as a
    single `Float64`, if there is only one mutational process. Multiple values indicate
    multiple simulataneous processes.
- `mutationdist::Vector{Symbol} = [:poisson]`: defines the distibution for new
    mutations (:poisson, :fixed, :poissontimedep, :fixedtimedep, :geometric). Length should
    match `length(μ)`.
- `ploidy::Int64 = 2`: cell ploidy
"""
struct MoranInput <: SinglelevelInput
    N::Int64
    tmax::Float64
    moranrate::Float64
    moranincludeself::Bool
    clonalmutations::Int64
    μ::Vector{Float64}
    mutationdist::Vector{Symbol}
    ploidy::Int64
end

function MoranInput(;
    N = 1000,
    tmax = 15.0,
    moranrate = 1.0,
    moranincludeself = true,
    clonalmutations = 0,
    μ = [1.0],
    mutationdist = fill(:poisson, length(μ)),
    ploidy = 2
)
    μ = tovector(μ)
    mutationdist = tovector(mutationdist)
    @assert length(μ) == length(mutationdist) "μ and mutationdist are not same length"

    return MoranInput(
        N, tmax, moranrate, moranincludeself, clonalmutations, μ, mutationdist, ploidy
    )
end

"""
    BranchingMoranInput <: SinglelevelInput

Input type for a single level simulation that grows by a branching process to `Nmax` cells
    and then switches to a Moran process.

# Fields:
- `N::Int64 = 1000`: number of cells
- `tmax::Float64 = 15.0`: maximum time to run simulation
- `moranrate::Float64 = 1.0`: moran update rate for wild-type cells
- `moranincludeself::Bool = true`: determines whether the same cell can be chosen to both
    divide and die in a moran step (in which case one offspring is killed)
- `birthrate::Float64 = moranrate`: birth rate for wild-type cells
- `deathrate::Float64 = 0.0`: death rate for wild-type cells
- `clonalmutations::Int64 = 0`: number of mutations shared by all cells
- `μ::Vector{Float64} = [1.0]`: mutation rate per division per cell. Can be passed as a
    single `Float64`, if there is only one mutational process. Multiple values indicate
    multiple simulataneous processes.
- `mutationdist::Vector{Symbol} = [:poisson]`: defines the distibution for new
    mutations (:poisson, :fixed, :poissontimedep, :fixedtimedep, :geometric). Length should
    match `length(μ)`.
- `ploidy::Int64 = 2`: cell ploidy
"""
struct BranchingMoranInput <: SinglelevelInput
    Nmax::Int64
    tmax::Float64
    moranrate::Float64
    moranincludeself::Bool
    birthrate::Float64
    deathrate::Float64
    clonalmutations::Int64
    μ::Vector{Float64}
    mutationdist::Vector{Symbol}
    ploidy::Int64
end

function BranchingMoranInput(;
    Nmax = 1000,
    tmax = 15.0,
    moranrate = 1.0,
    moranincludeself = true,
    birthrate = moranrate,
    deathrate = 0.0,
    clonalmutations = 0,
    μ = [1.0],
    mutationdist = fill(:poisson, length(μ)),
    ploidy = 2
)
    μ = tovector(μ)
    mutationdist = tovector(mutationdist)
    @assert length(μ) == length(mutationdist) "μ and mutationdist are not same length"

    return BranchingMoranInput(
        Nmax, tmax, moranrate, moranincludeself, birthrate, deathrate, clonalmutations, μ,
            mutationdist, ploidy
    )
end

"""
    newinput(::Type{InputType}; kwargs) where InputType <: SimulationInput

Create a new input of type `InputType`.
"""
function newinput(::Type{InputType}; kwargs) where InputType <: SimulationInput
    return InputType(;kwargs...)
end

"""
    newinput(input::InputType; kwargs...) where InputType <: SimulationInput

Create a new input of type `InputType`. Any fields not given in `kwargs` default to the
values in `input`.
"""
function newinput(input::InputType; kwargs...) where InputType <: SimulationInput
    newkwargs = Dict(
        field in keys(kwargs) ? field => kwargs[field] : field => getfield(input, field)
            for field in fieldnames(InputType))
    return InputType(;newkwargs...)
end


function Base.show(io::IO, input::BranchingInput)
    @printf(io, "Single level branching process:\n")
    @printf(io, "    Maximum cells = %d\n", input.Nmax)
    @printf(io, "    Maximum time = %.2f\n", input.tmax)
    @printf(
        io,
        "    Birth rate = %.3f, death rate = %.3f\n",
        input.birthrate,
        input.deathrate
    )
    for i in 1:length(input.μ)
        @printf(
            io,
            "    %s: μ = %.3f\n",
            mutationdist_string(input.mutationdist[i]),
            input.μ[i]
        )
    end
    @printf(io, "    Clonal mutations = %d\n", input.clonalmutations)
    @printf(io, "    Ploidy = %d\n", input.ploidy)
end

function Base.show(io::IO, input::MoranInput)
    @printf(io, "Single level Moran process:\n")
    @printf(io, "    Maximum cells = %d\n", input.N)
    @printf(io, "    Maximum time = %.2f\n", input.tmax)
    @printf(io, "    Moran rate = %.3f \n", input.moranrate)
    for i in 1:length(input.μ)
        @printf(
            io,
            "    %s: μ = %.3f\n",
            mutationdist_string(input.mutationdist[i]),
            input.μ[i]
        )
    end
    @printf(io, "    Clonal mutations = %d\n", input.clonalmutations)
    @printf(io, "    Ploidy = %d\n", input.ploidy)

end

function Base.show(io::IO, input::BranchingMoranInput)
    @printf(io, "Single level Branching -> Moran process:\n")
    @printf(io, "    Maximum cells = %d\n", input.Nmax)
    @printf(io, "    Maximum time = %.2f\n", input.tmax)
    @printf(
        io,
        "    Birth rate = %.3f, death rate = %.3f\n",
        input.birthrate,
        input.deathrate
    )
    @printf(io, "    Moran rate = %.3f\n", input.moranrate)
    for i in 1:length(input.μ)
        @printf(
            io, "    %s: μ = %.3f\n",
            mutationdist_string(input.mutationdist[i]),
            input.μ[i]
        )
    end
    @printf(io, "    Clonal mutations = %d\n", input.clonalmutations)
    @printf(io, "    Ploidy = %d\n", input.ploidy)

end

function mutationdist_string(mutationdist)
    if mutationdist == :fixed
        return "Fixed mutations"
    elseif mutationdist == :fixedtimedep
        return "Time-dependent fixed mutations"
    elseif mutationdist == :poisson
        return "Poisson distributed mutations"
    elseif mutationdist == :poissontimedep
        return "Time-dependent Poisson distibuted mutations"
    elseif mutationdist == :geometric
        return "Geometric distributed mutations"
    end
end
