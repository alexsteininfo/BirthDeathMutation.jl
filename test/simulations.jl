@testset "neutral simulate! BirthDeathBlock" begin
    rng = MersenneTwister(100)
    Nmax = 10
    block = BirthDeathBlock(
        birthrate    = (s, N) -> 1.0 * (1 + s),
        deathrate    = (s, N) -> 0.0,
        stopfunction = pop -> popsize(pop) >= Nmax,
        μ = 1.0,
        mutationdist = :poisson,
    )
    pop = initialize_population(1)
    pop = simulate!(pop, block, rng)
    @test length(allcells(pop)) == Nmax
end

@testset "neutral simulate! MoranBlock" begin
    rng = MersenneTwister(100)
    tmax = 10.0
    N = 10
    moran = MoranBlock(
        N            = N,
        moranrate    = (s, N) -> 1.0 * (1 + s),
        stopfunction = pop -> age(pop) >= tmax,
        μ = 1.0,
        mutationdist = :poisson,
    )
    pop = initialize_population(N)
    pop = simulate!(pop, moran, rng)
    @test length(allcells(pop)) == N
    @test age(pop) <= tmax + 1.0  # may slightly exceed tmax in last step
end

@testset "chained blocks" begin
    rng = MersenneTwister(42)
    Nmax = 20
    tmax = 5.0

    b1 = BirthDeathBlock(
        birthrate    = (s, N) -> 1.0 * (1 + s),
        deathrate    = (s, N) -> 0.0,
        stopfunction = pop -> popsize(pop) >= Nmax,
        μ = 1.0,
        mutationdist = :poisson,
    )
    b2 = MoranBlock(
        N            = Nmax,
        moranrate    = (s, N) -> 1.0 * (1 + s),
        stopfunction = pop -> age(pop) >= tmax,
        μ = 1.0,
        mutationdist = :poisson,
    )

    pop = initialize_population(1)
    t_after_b1 = let pop = simulate!(pop, b1, rng); age(pop) end
    pop = initialize_population(1)
    pop = simulate!(pop, b1, rng)
    @test length(allcells(pop)) == Nmax
    t_mid = age(pop)

    pop = simulate!(pop, b2, rng)
    @test length(allcells(pop)) == Nmax
    @test age(pop) >= t_mid  # time is monotonically increasing
end

@testset "selection simulate! BirthDeathBlock" begin
    rng = MersenneTwister(100)
    Nmax = 20
    block = BirthDeathBlock(
        birthrate    = (s, N) -> 1.0 * (1 + s),
        deathrate    = (s, N) -> 0.0,
        stopfunction = pop -> popsize(pop) >= Nmax,
        selection    = SelectionPredefined([0.5], [0.5]),
        μ = 1.0,
        mutationdist = :poisson,
    )
    pop = initialize_population(1)
    pop = simulate!(pop, block, rng)
    @test length(allcells(pop)) == Nmax
    @test sum(getsubclonesizes(pop)) == Nmax
    @test getsubclonesizes(pop) == counts(getclonetype.(allcells(pop)), 1:length(pop.subclones))
end

@testset "clone-spec initialization" begin
    pop = initialize_population([
        (n=9, s=0.0, mutations=50),
        (n=1, s=0.1, mutations=55),
    ])
    @test length(allcells(pop)) == 10
    @test length(pop.subclones) == 2
    @test pop.subclones[1].s == 0.0
    @test pop.subclones[2].s == 0.1
    @test pop.subclones[1].size == 9
    @test pop.subclones[2].size == 1
    wt_muts = [allcells(pop)[i].data.mutations for i in 1:9]
    @test all(wt_muts .== 50)
    fit_cell = allcells(pop)[10]
    @test fit_cell.data.mutations == 55
    @test fit_cell.data.clonetype == 2
end

@testset "population size mismatch error" begin
    rng = MersenneTwister(1)
    moran = MoranBlock(
        N            = 5,
        moranrate    = (s, N) -> 1.0,
        stopfunction = pop -> age(pop) >= 1.0,
    )
    pop = initialize_population(10)  # wrong size
    @test_throws ErrorException simulate!(pop, moran, rng)
end

@testset "restart_on_extinction" begin
    # Default is false
    block_default = BirthDeathBlock(
        birthrate    = (s, N) -> 1.0,
        deathrate    = (s, N) -> 0.5,
        stopfunction = pop -> popsize(pop) >= 5,
    )
    @test block_default.restart_on_extinction == false

    # With restart=true and net-positive dynamics, always reaches target
    rng = MersenneTwister(99)
    block = BirthDeathBlock(
        birthrate             = (s, N) -> 1.0,
        deathrate             = (s, N) -> 0.5,
        stopfunction          = pop -> popsize(pop) >= 10,
        restart_on_extinction = true,
    )
    pop = initialize_population(1)
    pop = simulate!(pop, block, rng)
    @test popsize(pop) >= 10
    @test age(pop) >= 0.0
end
