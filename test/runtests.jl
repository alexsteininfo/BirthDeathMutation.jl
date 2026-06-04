using BirthDeathMutation
using Test
using Random
using StatsBase
using AbstractTrees
using Distributions

tests = [
    "initialisation",
    "simulations",
    "process_mutations",
    "treesimulations",
]

@testset "BirthDeathMutation.jl" begin
    for test in tests
        @testset "$test" begin
            include(test*".jl")
        end
    end
end
