# BirthDeathMutation

Julia package for simulating single-population somatic evolution with full lineage tracking.

This project started from [SomaticEvolution.jl](https://github.com/jessierenton/SomaticEvolution.jl) but has been substantially simplified: it supports only `SimpleTreeCell` (tree-based lineage tracking), uses a composable **block** API instead of monolithic input types, and removes multilevel dynamics, VAF analysis, and time-series wrappers.

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/alexanderstein/BirthDeathMutation.jl")
```

## How it works

Simulations are built from two composable block types:

| Block | Description |
|---|---|
| `BirthDeathBlock` | Stochastic branching process; birth and death rates are arbitrary functions of `s` and `N` |
| `MoranBlock` | Moran process at fixed population size `N`; Moran rate is a function of `s` and `N` |

Each block runs until its `stopfunction(pop) -> Bool` returns `true`. Blocks are applied
imperatively with `simulate!`, which modifies the population in-place and returns it:

```julia
pop = initialize_population(1)         # one founding cell
pop = simulate!(pop, block1, rng)      # extend tree with block1
pop = simulate!(pop, block2, rng)      # extend tree with block2 (starts from block1's result)
```

### Block parameters

```julia
BirthDeathBlock(;
    birthrate,       # b(s::Float64, N::Int) -> Float64
    deathrate,       # d(s::Float64, N::Int) -> Float64
    stopfunction,    # stop(pop::Population) -> Bool
    selection  = NeutralSelection(),
    μ          = 1.0,
    mutationdist = :poisson,      # :poisson | :fixed | :geometric
    restart_on_extinction = false # if true, retry from initial state on extinction
)

MoranBlock(;
    N,               # fixed population size (must match incoming pop)
    moranrate,       # r(s::Float64, N::Int) -> Float64
    stopfunction,    # stop(pop::Population) -> Bool
    selection  = NeutralSelection(),
    μ          = 1.0,
    mutationdist = :poisson,
    moranincludeself = false
)
```

### Selection

Selection is specified per block via the `selection` keyword:

| Type | Description |
|---|---|
| `NeutralSelection()` | All cells identical fitness (default) |
| `SelectionPredefined(coefficients, times)` | Fit mutants arise at specified times with specified `s` values |
| `SelectionDistribution(dist, probability, max_subclones)` | Fit mutants arise per division with `s` drawn from `dist` |

The selection coefficient `s` is stored on each `Subclone` and passed directly to the
block's rate functions, so the user controls how fitness enters the dynamics.

### Initial conditions

```julia
# N identical cells, no initial mutations
pop = initialize_population(N)

# N cells with n_mutations clonal mutations
pop = initialize_population(N; n_mutations=50)

# Multiple clones with different s and mutation counts
pop = initialize_population([
    (n=900, s=0.0, mutations=50),
    (n=100, s=0.1, mutations=55),
])

# Output of a prior simulate! call is directly usable as input
pop = simulate!(pop, block1, rng)  # pop is now the initial condition for block2
pop = simulate!(pop, block2, rng)
```

## Cell representation

Cells are stored as `BinaryNode{SimpleTreeCell}` nodes. Full ancestry is preserved
throughout the simulation as a binary tree, so any pair of cells can be traced to their
MRCA. Each cell stores the *number* of mutations acquired since its parent divided.

## Examples

Three worked examples are in the `examples/` directory:

| File | Description |
|---|---|
| `SingleCellExpansion.jl` | Logistic growth from one cell to K=10 000 with stochastic selection; uses `restart_on_extinction = true` |
| `ChainedBlocks.jl` | Three chained phases: exponential growth → Moran homeostasis → logistic re-expansion |
| `ArbitraryInitialCondition.jl` | Declining population (b < d) starting from 5 pre-existing clones with distinct fitness values |

Run any example with:

```julia
julia --project=. examples/SingleCellExpansion.jl
```

Minimal inline example — logistic growth from one cell with automatic restart on extinction:

```julia
using BirthDeathMutation, Random, Distributions

rng = MersenneTwister(42)
K = 10_000

block = BirthDeathBlock(
    birthrate             = (s, N) -> 1.0 * (1 + s),
    deathrate             = (s, N) -> 0.8 * N / K,
    stopfunction          = pop -> popsize(pop) >= K,
    selection             = SelectionDistribution(Exponential(0.1), 0.01, 3),
    μ                     = 1.0,
    restart_on_extinction = true,   # retry if lineage goes extinct
)

pop = initialize_population(1)
pop = simulate!(pop, block, rng)
println("cells: ", popsize(pop))
println("mean mutations: ", average_mutations(pop))
```

## Analysis functions

- `mutations_per_cell(pop)` — mutations per alive cell (summed along lineage)
- `average_mutations(pop)` — mean mutations across all alive cells
- `clonal_mutations(pop)` — mutations shared by every alive cell (at MRCA)
- `pairwise_differences(pop[, idx])` — pairwise distance frequency dict
- `pairwisedistances(pop[, idx])` — raw vector of pairwise distances
- `pairwisedistance(cell1, cell2)` — distance between two specific cells
- `coalescence_times(pop[, idx])` — time-to-MRCA for every pair of alive cells
- `getsubclonesizes(pop)` — vector of cell counts per subclone

## Customisation

New block types can be added by defining a new struct and a corresponding `simulate!`
method dispatching on it. New selection regimes can be added by defining new
`AbstractSelection` subtypes with `newsubclone_ready` and `getselectioncoefficient`
methods. See `src/blocks.jl` and `src/selection.jl` for reference.
