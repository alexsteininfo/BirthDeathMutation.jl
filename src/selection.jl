"""
    AbstractSelection

Type that defines how and when fit mutants will arise in the simulation.

# Subtypes
- `NeutralSelection`
- `SelectionPredefined`
- `SelectionDistribution`
"""
abstract type AbstractSelection end

"""
    NeutralSelection <: AbstractSelection

Neutral selection type: all mutations have the same fitness.
"""
struct NeutralSelection <: AbstractSelection
end

"""
    SelectionPredefined <: AbstractSelection

Non-neutral selection type: fit mutations arise at given times and with given selection 
coefficients.

# Fields
- `mutant_selection::Vector{Float64}` -- selection coefficient for each fit mutant that 
    arises
- `mutant_selection::Vector{Int64}` -- time that each fit mutant arises
"""
struct SelectionPredefined <: AbstractSelection
    mutant_selection::Vector{Float64}
    mutant_time::Vector{Float64}
end

"""
    SelectionDistribution{D<:Distribution} <: AbstractSelection

Non-neutral selection type: fit mutations arise at each cell division with a given 
    probability and with selection coefficient drawn from the distribution, up to a maximum 
    number.

# Fields
- `distribution::D` -- used to draw selection coefficients using `rand([rng, ]distibution)`
- `mutation_probability::Float64` -- probability that a fit mutant arises during cell 
    division
- `maximum_subclones::Int64` -- maximum number of fit mutants that can occur

"""
struct SelectionDistribution{D<:Distribution} <:AbstractSelection
    distribution::D
    mutation_probability::Float64
    maximum_subclones::Int64
end

"""
    getmaxsubclones(selection::AbstractSelection)

Returns the maximum number of subclones that can be generated from `selection`.
"""
function getmaxsubclones end

getmaxsubclones(::NeutralSelection) = 1
getmaxsubclones(selection::SelectionPredefined) = length(selection.mutant_selection) + 1
getmaxsubclones(selection::SelectionDistribution) = selection.maximum_subclones

"""
    newsubclone_ready(selection::AbstractSelection, nsubclonescurrent, nsubclones, t, rng)

Returns `true` if a new fit mutant (subclone) should be introduced.
"""
function newsubclone_ready end

newsubclone_ready(::NeutralSelection, nsubclonescurrent, nsubclones, t, rng) = false

function newsubclone_ready(selection::SelectionPredefined, nsubclonescurrent, nsubclones, t, rng)
    return nsubclonescurrent < nsubclones && t >= selection.mutant_time[nsubclonescurrent]
end

function newsubclone_ready(selection::SelectionDistribution, nsubclonescurrent, nsubclones, t, rng)
    return nsubclonescurrent < nsubclones &&  rand(rng) < selection.mutation_probability 
end

"""
    SelectionAccumulative{D<:Distribution} <: AbstractSelection

Multiplicative fitness accumulation with diminishing-returns epistasis (eq 17,
Leventhal et al. Nature Communications 2026).

Each new subclone's birth rate is computed from its parent's birth rate `r_p`:

    r_new = min( r_p * (1 + s·X·(1 - r_p/M)),  M )

where `X` is drawn from `factor_distribution` (use `Dirac(1.0)` for fixed effects,
`Exponential(1.0)` for random effects with unit-mean factor).

# Fields
- `s::Float64` — fitness step parameter
- `M::Float64` — maximum birth rate cap
- `factor_distribution::D` — distribution of the random factor X
- `mutation_probability::Float64` — per-division probability of a fit mutation
- `maximum_subclones::Int64` — cap on the number of tracked subclones
"""
struct SelectionAccumulative{D<:Distribution} <: AbstractSelection
    s::Float64
    M::Float64
    factor_distribution::D
    mutation_probability::Float64
    maximum_subclones::Int64
end

getmaxsubclones(sel::SelectionAccumulative) = sel.maximum_subclones

function newsubclone_ready(sel::SelectionAccumulative, nsubclonescurrent, nsubclones, _, rng)
    return nsubclonescurrent < nsubclones && rand(rng) < sel.mutation_probability
end

"""
    SelectionMaxRandom{D<:Distribution} <: AbstractSelection

Random fitness effects, non-multiplicative. A new subclone's selection coefficient
is `max(s_parent, draw)` where `draw ~ distribution`.  This means a child subclone
is never less fit than its parent.

# Fields
- `distribution::D` — distribution from which new fitness effects are drawn
- `mutation_probability::Float64` — per-division probability of a fit mutation
- `maximum_subclones::Int64` — cap on the number of tracked subclones
"""
struct SelectionMaxRandom{D<:Distribution} <: AbstractSelection
    distribution::D
    mutation_probability::Float64
    maximum_subclones::Int64
end

getmaxsubclones(sel::SelectionMaxRandom) = sel.maximum_subclones

function newsubclone_ready(sel::SelectionMaxRandom, nsubclonescurrent, nsubclones, t, rng)
    return nsubclonescurrent < nsubclones && rand(rng) < sel.mutation_probability
end

"""
    getselectioncoefficient(selection::AbstractSelection, nsubclonescurrent, rng)

Returns the next selection coefficient for a new mutant.
The optional keyword `parent_s` is the selection coefficient of the parent subclone;
it is used by types that accumulate fitness across generations.
"""
function getselectioncoefficient end

function getselectioncoefficient(selection::SelectionPredefined, nsubclonescurrent, rng)
    selection.mutant_selection[nsubclonescurrent]
end
function getselectioncoefficient(selection::SelectionDistribution, nsubclonescurrent, rng)
    return rand(rng, selection.distribution)
end

# Backwards-compatible fallback: existing types ignore parent_s
function getselectioncoefficient(sel::AbstractSelection, n, rng; parent_s=0.0)
    getselectioncoefficient(sel, n, rng)
end

function getselectioncoefficient(sel::SelectionAccumulative, _n, rng; parent_s=1.0)
    X     = rand(rng, sel.factor_distribution)
    r_new = parent_s * (1 + sel.s * X * (1 - parent_s / sel.M))
    return min(r_new, sel.M)
end

function getselectioncoefficient(sel::SelectionMaxRandom, _n, rng; parent_s=0.0)
    draw = rand(rng, sel.distribution)
    return max(parent_s, draw)
end
