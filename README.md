# BirthDeathMutation

Julia package for simulating single-population somatic evolution with full lineage tracking.

This project started from [SomaticEvolution.jl](https://github.com/jessierenton/SomaticEvolution.jl) but has been substantially modified: multilevel (multi-module) dynamics have been removed, only `SimpleTreeCell` (tree-based lineage tracking) is supported, and the package has been renamed to reflect a single birth-death-mutation population focus.

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/alexanderstein/BirthDeathMutation")
```

## How it works

Simulations are run using `runsimulation`, `runsimulation_timeseries`, or
`runsimulation_timeseries_returnfinalpop`.

First create an input that defines the population dynamics:

| Input type | Description |
|---|---|
| `BranchingInput` | Branching process from a single cell until `Nmax` cells |
| `MoranInput` | Moran process starting from `N` identical cells |
| `BranchingMoranInput` | Branching process to `Nmax`, then switches to Moran |

All input types share common parameters: `μ` (mutation rate per division), `mutationdist`
(`:poisson`, `:fixed`, `:poissontimedep`, `:fixedtimedep`, `:geometric`),
`clonalmutations`, and `ploidy`.

By default selection is neutral. Other regimes can be specified with an `AbstractSelection`
argument:
- `NeutralSelection()` — default, all cells identical fitness
- `SelectionPredefined(mutant_selection, mutant_time)` — fit mutants at specified times
- `SelectionDistribution(dist, probability, max_subclones)` — fit mutants arise
  stochastically with selection coefficients drawn from `dist`

## Cell representation

Cells are stored as `BinaryNode{SimpleTreeCell}` nodes. Full ancestry is preserved
throughout the simulation as a binary tree, so any pair of cells can be traced to their
MRCA. Each cell stores the *number* of mutations acquired since its parent divided (not
individual mutation IDs).

## Examples

Branching process from a single cell until 100 cells:
```julia
using BirthDeathMutation, Random

input = BranchingInput(Nmax=100)
rng = Random.seed!(12)
simulation = runsimulation(input, rng)
```

Moran process with selection:
```julia
using BirthDeathMutation, Random, Distributions

input = MoranInput(N=500, tmax=20.0)
selection = SelectionDistribution(Exponential(0.2), 0.1, 10)
rng = Random.seed!(42)
simulation = runsimulation(input, selection, rng)
```

Branching process recording population state at multiple timepoints:
```julia
using BirthDeathMutation, Random

input = BranchingInput(Nmax=1000)
rng = Random.seed!(1)
timesteps = 1.0:1.0:10.0
data = runsimulation_timeseries(input, timesteps, average_mutations, rng)
```

## Analysis functions

- `mutations_per_cell(simulation)` — mutations per alive cell (summed along lineage)
- `average_mutations(simulation)` — mean mutations across all alive cells
- `clonal_mutations(simulation)` — mutations shared by every alive cell (at MRCA)
- `pairwise_differences(simulation)` — pairwise distance frequency dict
- `pairwisedistances(simulation)` — raw vector of pairwise distances
- `coalescence_times(root)` — time-to-MRCA for every pair of alive cells
- `time_to_MRCA(node1, node2, t)` — time since MRCA of two specific cells

## Customisation

New simulation dynamics can be added by defining new `SimulationInput` subtypes and
corresponding `simulate!` and `initialize_population` methods. New selection regimes can
be added by defining new `AbstractSelection` subtypes. See the existing implementations
in `src/` for examples.
