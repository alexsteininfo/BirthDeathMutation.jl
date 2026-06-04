@testset "mutation distributions" begin
    rng = MersenneTwister(42)

    # :fixed always returns the rounded mean
    @test BirthDeathMutation.numbernewmutations(rng, :fixed, 3.0) == 3

    # :poisson mean ≈ μ over many draws
    samples = [BirthDeathMutation.numbernewmutations(rng, :poisson, 2.0) for _ in 1:1000]
    @test abs(mean(samples) - 2.0) < 0.2

    # :geometric mean ≈ μ over many draws
    samples = [BirthDeathMutation.numbernewmutations(rng, :geometric, 3.0) for _ in 1:1000]
    @test abs(mean(samples) - 3.0) < 0.5

    # Unknown distribution throws
    @test_throws ErrorException BirthDeathMutation.numbernewmutations(rng, :unknown, 1.0)
end
