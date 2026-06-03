"""
    BirthDeathBlock{F1,F2,F3,S}

Defines a birth-death (branching) simulation block. Birth and death rates are
user-provided functions of selection coefficient `s` and population size `N`.
The block runs until `stopfunction(pop)` returns `true`.

# Fields
- `birthrate::F1` — `b(s::Float64, N::Int) -> Float64`
- `deathrate::F2` — `d(s::Float64, N::Int) -> Float64`
- `stopfunction::F3` — `stop(pop::Population) -> Bool`
- `selection::S` — how fit mutants arise (default: `NeutralSelection()`)
- `μ::Float64` — mean mutations per division (default: `1.0`)
- `mutationdist::Symbol` — mutation distribution: `:poisson`, `:fixed`, or `:geometric`
- `restart_on_extinction::Bool` — if `true`, restart from the initial population
  whenever the population goes extinct (default: `false`)
"""
struct BirthDeathBlock{F1<:Function, F2<:Function, F3<:Function, S<:AbstractSelection}
    birthrate::F1
    deathrate::F2
    stopfunction::F3
    selection::S
    μ::Float64
    mutationdist::Symbol
    restart_on_extinction::Bool
end

function BirthDeathBlock(;
    birthrate,
    deathrate,
    stopfunction,
    selection             = NeutralSelection(),
    μ                     = 1.0,
    mutationdist          = :poisson,
    restart_on_extinction = false
)
    return BirthDeathBlock(
        birthrate, deathrate, stopfunction,
        selection, Float64(μ), mutationdist, restart_on_extinction
    )
end

"""
    MoranBlock{F1,F2,S}

Defines a Moran process simulation block with fixed population size `N`.
The Moran rate is a user-provided function of selection coefficient `s` and
population size `N`. The block runs until `stopfunction(pop)` returns `true`.

# Fields
- `N::Int` — fixed population size (population passed to `simulate!` must match)
- `moranrate::F1` — `r(s::Float64, N::Int) -> Float64`
- `stopfunction::F2` — `stop(pop::Population) -> Bool`
- `selection::S` — how fit mutants arise (default: `NeutralSelection()`)
- `μ::Float64` — mean mutations per division (default: `1.0`)
- `mutationdist::Symbol` — mutation distribution: `:poisson`, `:fixed`, or `:geometric`
- `moranincludeself::Bool` — whether a cell can replace itself (default: `false`)
"""
struct MoranBlock{F1<:Function, F2<:Function, S<:AbstractSelection}
    N::Int
    moranrate::F1
    stopfunction::F2
    selection::S
    μ::Float64
    mutationdist::Symbol
    moranincludeself::Bool
end

function MoranBlock(;
    N,
    moranrate,
    stopfunction,
    selection = NeutralSelection(),
    μ = 1.0,
    mutationdist = :poisson,
    moranincludeself = false
)
    return MoranBlock(N, moranrate, stopfunction, selection, Float64(μ), mutationdist, moranincludeself)
end
