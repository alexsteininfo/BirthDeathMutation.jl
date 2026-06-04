##
## Example 1: Single-cell logistic expansion with stochastic selection
##
## A single founding cell grows logistically to carrying capacity K = 10 000.
## The death rate rises linearly with population size (d(N) = d0·N/K), so net
## growth decelerates as N → K. Fit mutants arise stochastically at each
## division (probability 0.01) with selection coefficients drawn from
## Exponential(0.1). Because a single founding cell has a high probability of
## stochastic extinction, restart_on_extinction = true retries automatically
## until the lineage survives.
##

using Pkg
Pkg.activate(dirname(@__DIR__))

using BirthDeathMutation
using Random
using Distributions
using Statistics: mean

rng = MersenneTwister(12)

K  = 10_000   # carrying capacity
b  = 1.0      # intrinsic birth rate
d0 = 0.8      # death rate at N = K (logistic equilibrium: b·(1+s) ≈ d0 when s=0)

block = BirthDeathBlock(
    birthrate             = (s, N) -> b * (1 + s),
    deathrate             = (s, N) -> d0 * N / K,
    stopfunction          = pop -> popsize(pop) >= K,
    selection             = SelectionDistribution(Exponential(0.1), 0.01, 3),
    μ                     = 1.0,
    mutationdist          = :poisson,
    restart_on_extinction = true,
)

pop = initialize_population(1)
pop = simulate!(pop, block, rng)

println("Population size:     ", popsize(pop))
println("Simulation time:     ", round(age(pop), digits=3))
println("Number of subclones: ", length(pop.subclones))
for (i, sc) in enumerate(pop.subclones)
    println("  Subclone $i — s = $(round(sc.s, digits=3)), size = $(sc.size)")
end

mpc = mutations_per_cell(pop)
println("\nMutations per cell:  mean=$(round(mean(mpc), digits=1)), " *
        "min=$(minimum(mpc)), max=$(maximum(mpc))")
println("Clonal mutations:    ", clonal_mutations(pop))

pd = pairwise_differences(pop, 1:min(50, popsize(pop)))
println("Pairwise differences (sample of 50 cells):")
for (d, n) in sort(collect(pd))
    println("  $d differences: $n pairs")
end
