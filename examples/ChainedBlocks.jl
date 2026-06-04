##
## Example 2: Chained blocks — growth, homeostasis, and a second expansion
##
## Demonstrates the composable block API:
##
##   Phase 1  BirthDeathBlock   Exponential growth from 1 cell to 500,
##                               with predefined selection event at t ≈ 0.5
##
##   Phase 2  MoranBlock        Moran homeostasis at N = 500 for 30 time units,
##                               neutral (no new subclones)
##
##   Phase 3  BirthDeathBlock   Second growth phase from 500 → 2000 cells
##                               with logistic-like death rate (d depends on N)
##
## Phase 2's stopfunction depends on the time at which phase 1 ends, so it is
## constructed after phase 1 completes. The tree is extended in-place across all
## three phases.
##

using Pkg
Pkg.activate(dirname(@__DIR__))

using BirthDeathMutation
using Random
using Statistics: mean

rng = MersenneTwister(7)

# --- Block definitions ---

phase1 = BirthDeathBlock(
    birthrate    = (s, N) -> 1.0 * (1 + s),
    deathrate    = (s, N) -> 0.0,
    stopfunction = pop -> popsize(pop) >= 500,
    selection    = SelectionPredefined([0.3], [0.5]),
    μ            = 1.0,
    mutationdist = :poisson,
)

K = 2000
phase3 = BirthDeathBlock(
    birthrate    = (s, N) -> 1.0 * (1 + s),
    deathrate    = (s, N) -> 0.8 * N / K,
    stopfunction = pop -> popsize(pop) >= K,
    μ            = 1.0,
    mutationdist = :poisson,
)

function print_status(label, pop)
    println(label)
    println("  N = $(popsize(pop)), t = $(round(age(pop), digits=2))")
    println("  subclones: $(length(pop.subclones)) — sizes: $(getsubclonesizes(pop))")
end

# --- Run ---

pop = initialize_population(1)

pop = simulate!(pop, phase1, rng)
t0 = age(pop)
print_status("After phase 1 (growth to 500):", pop)

# Phase 2 is defined here so its stopfunction can capture t0
phase2 = MoranBlock(
    N            = 500,
    moranrate    = (s, N) -> 1.0 * (1 + s),
    stopfunction = let t = t0; pop -> age(pop) >= t + 30.0; end,
    μ            = 0.5,
    mutationdist = :poisson,
)

pop = simulate!(pop, phase2, rng)
print_status("After phase 2 (Moran 30 t.u.):", pop)

pop = simulate!(pop, phase3, rng)
print_status("After phase 3 (growth to 2000):", pop)

# --- Analysis ---
println("\nAnalysis of final population:")
mpc = mutations_per_cell(pop)
println("  Mean mutations per cell: ", round(mean(mpc), digits=1))
println("  Clonal mutations:        ", clonal_mutations(pop))

idx = randperm(rng, popsize(pop))[1:min(30, popsize(pop))]
ct = coalescence_times(pop, idx)
println("  Mean coalescence time (sample of 30): ", round(mean(ct), digits=2))
