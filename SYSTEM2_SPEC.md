# System 2: High-Rate Driver Mutation Package — Specification

## Motivation

`BirthDeathMutation.jl` (System 1) tracks individual fit clones as named subclones in a cell
lineage tree. This works well when the driver mutation probability ν is low (0.001–0.5), so that
clones are sparse, identifiable, and competitively meaningful. At high ν (≥ 1), every cell
acquires a driver every few divisions — clone identities lose meaning and the subclone tree
explodes in size.

System 2 is a **separate Julia package** aimed at the high-ν regime. It focuses on:
- The distribution of accumulated driver mutations per cell over time
- How selection shapes the fitness landscape of the population
- Dynamics at the level of mutation-count classes, not individual clones

---

## Core model

### State

A population of N cells. Each cell carries an integer count of driver mutations, `k ∈ {0, 1, 2, ...}`.

There are **no neutral mutations** (μ = 0) and **no clone identities** — only the fitness
landscape n[k] (number of cells with exactly k drivers) matters.

### Fitness

Each cell's birth rate and death rate are functions of its driver count k and total population N:

```
birthrate(k, N)    # user-supplied; e.g. (k, N) -> 1 + k*s
deathrate(k, N)    # user-supplied; e.g. (k, N) -> 0.0
```

### Driver mutation model

Each daughter cell independently acquires j ~ Poisson(ν) driver mutations per division.
Each driver draws a fitness increment from a distribution D (e.g., Exponential(s)).

Because fitness increments are **random per mutation event**, two cells with the same k can
have different accumulated fitnesses — the simple bin-by-k structure breaks down once fitness is
random.

---

## Per-cell tracking: the required architecture

Since fitness increments are random, cell fitness must be tracked individually.

### Cell representation

```julia
struct DriverCell
    k::Int      # accumulated driver mutation count
    fitness::Float64   # accumulated fitness (birthrate factor)
end
```

Population state: `Vector{DriverCell}` of length N.

### Gillespie loop

At each step:
1. Compute `Rmax = maximum(birthrate(c.fitness, N) + deathrate(c.fitness, N) for c in cells)`
   — O(N) per step.
2. Sample waiting time `Δt ~ Exp(1) / (Rmax × N)`.
3. Pick a random cell.
4. Accept/reject birth or death via the usual thinning step.
5. On birth: create two daughter cells. For each daughter, draw j ~ Poisson(ν); for each of
   the j mutations, draw δ from distribution D; update fitness accordingly (additive, multiplicative,
   or via a user-supplied update rule).
6. On death: remove the cell.

**Complexity**: O(N) per step → O(N²) total for a simulation reaching size N from 1 cell.
For N = 10,000 this is ~10⁸–10⁹ operations in Julia — feasible (seconds to minutes), but not fast.

### Optimisation options

- **Maintain a running max**: track `current_rmax` and update it incrementally on birth/death.
  Only recompute from scratch when a cell with `fitness == current_rmax` dies. This reduces
  Rmax computation to O(1) most steps.
- **Sorted structure**: keep cells sorted by fitness (a heap). O(log N) update per step.
- **Vectorised rejection**: batch multiple steps using tau-leaping for approximate but fast
  large-N behaviour.

---

## Suggested API

```julia
using DriverLoadSim   # working name for the package

result = simulate(;
    N_target      = 10_000,
    ν             = 2.0,           # Poisson driver rate per daughter cell
    fitness_init  = 1.0,           # starting fitness of the founding cell
    fitness_update = (f, δ) -> f + δ,   # how each driver changes fitness
    driver_dist   = Exponential(0.1),   # distribution of fitness increments
    birthrate     = (f, N) -> f,        # birth rate as function of fitness and N
    deathrate     = (f, N) -> 0.0,
    trajectory_dt = 0.5,           # time spacing for trajectory recording
    rng           = Random.GLOBAL_RNG,
)
```

### Output

```julia
struct DriverLoadResult
    trajectory::Vector{DriverLoadPoint}   # time series
    final_cells::Vector{DriverCell}       # state at end
end

struct DriverLoadPoint
    t::Float64
    N_total::Int
    mean_k::Float64          # mean driver count
    var_k::Float64           # variance in driver count
    mean_fitness::Float64
    var_fitness::Float64
    k_distribution::Vector{Int}   # histogram: k_distribution[i] = cells with k=i-1
end
```

---

## Analysis focus

Unlike System 1 (which focuses on SFS, clone sizes, and neutral mutations per cell), System 2
analysis will focus on:

- **Fitness distribution over time**: how does the fitness histogram change?
- **Mean and variance of fitness**: does variance track the mutation-selection balance?
- **Driver load distribution**: distribution of k values, analogous to the neutral SFS
- **Selection signatures**: does the tail of the fitness distribution grow faster than expected
  under neutrality?
- **Mutation accumulation rate**: does mean k grow as ν × 2t (neutral) or faster (positive selection)?

---

## Relationship to System 1

| Feature                   | System 1 (`BirthDeathMutation.jl`) | System 2 (new package) |
|---------------------------|------------------------------------|------------------------|
| ν range                   | 0.001–0.5 (probability mode)       | 0.1–4.0+               |
| Clone tracking            | Yes (subclone tree)                | No                     |
| Neutral mutations (μ)     | Yes                                | No                     |
| Fitness per cell          | Per subclone (shared)              | Per cell (individual)  |
| State representation      | Cell tree + subclones              | Vector of DriverCells  |
| Gillespie complexity      | O(|subclones|) per step            | O(N) per step          |
| Analysis outputs          | SFS, clone sizes, mut/cell         | Fitness distribution, k histogram |

The two packages are complementary — System 1 for studying clonal dynamics with rare drivers,
System 2 for studying fitness landscape evolution under frequent drivers.

---

## Implementation notes

- Start from scratch: no shared code with `BirthDeathMutation.jl` beyond the Gillespie pattern.
- The `fitness_update` function replaces the selection type hierarchy of System 1.
- No `Subclone`, `SimpleTreeCell`, or `BinaryNode` needed — `DriverCell` is a flat struct.
- `restart_on_extinction` logic is identical to System 1 and can be copied.
- Tests: check that at ν=0, mean fitness stays constant; at ν>0 with positive selection, mean
  fitness increases over time; at ν=0.001, results agree approximately with System 1.
