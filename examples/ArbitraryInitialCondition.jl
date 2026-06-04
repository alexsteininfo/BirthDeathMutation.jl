##
## Example 3: Declining population from an arbitrary initial condition
##
## Starts from 10 000 cells distributed across 5 pre-existing clones. Each clone
## has a distinct selection coefficient; mutations per cell are drawn independently
## from a discrete uniform distribution U{0,…,10}. Because birthrate < deathrate
## the total population declines. The simulation stops the first time the
## population reaches N_min = 1 000. Fitter clones decline more slowly and
## therefore represent a larger fraction of the surviving cells.
##

using Pkg
Pkg.activate(dirname(@__DIR__))

using BirthDeathMutation
using Random
using Statistics: mean

rng = MersenneTwister(42)

# --- Initial condition: 5 clones, equal size, distinct fitness ---

N_per_clone = 2_000
N_min       = 1_000
clone_s     = [-0.1, -0.05, 0.0, 0.05, 0.1]

# Mutations per cell: one independent draw from U{0, ..., 10} per clone
clone_mutations = rand(rng, 0:10, 5)

clones = [(n=N_per_clone, s=clone_s[i], mutations=clone_mutations[i]) for i in 1:5]

println("Initial clone specification (total N = $(5 * N_per_clone)):")
for (i, c) in enumerate(clones)
    println("  Clone $i — s=$(c.s), cells=$(c.n), mutations/cell=$(c.mutations)")
end

pop = initialize_population(clones)

# --- Simulation: net-negative growth until first passage to N_min ---

block = BirthDeathBlock(
    birthrate    = (s, N) -> 0.5 * (1 + s),
    deathrate    = (s, N) -> 0.8,
    stopfunction = pop -> popsize(pop) <= N_min,
    μ            = 1.0,
    mutationdist = :poisson,
)

simulate!(pop, block, rng)

# --- Results ---

println("\nAfter first passage to N_min = $N_min:")
println("  Final population size: ", popsize(pop))
println("  Simulation time:       ", round(age(pop), digits=3))
sizes = getsubclonesizes(pop)
println("  Clone composition:")
for (i, sc) in enumerate(pop.subclones)
    frac = round(100 * sizes[i] / popsize(pop), digits=1)
    println("  Clone $i — s=$(sc.s), cells=$(sizes[i]) ($frac%)")
end

mpc = mutations_per_cell(pop)
println("\nMutations per cell: mean=$(round(mean(mpc), digits=1)), " *
        "min=$(minimum(mpc)), max=$(maximum(mpc))")
