@testset "neutral runsimulation" begin
    rng = MersenneTwister(100)
    input = BranchingInput(
        Nmax=10,
        mutationdist=:poisson,
        birthrate=1,
        deathrate=0.0,
        clonalmutations=0,
        μ=1
    )
    simulation = runsimulation(input, rng)
    @test length(allcells(simulation.output)) == 10

    tmax = 10
    input = MoranInput(
        N=10,
        tmax=tmax,
        mutationdist=:poisson,
        moranrate=1.0,
        clonalmutations=0,
        μ=1
    )
    simulation = runsimulation(input, rng)
    @test length(allcells(simulation.output)) == 10
    @test age(simulation) == simulation.output.t
    @test age(simulation) <= tmax

    tmax = 10
    input = BranchingMoranInput(
        Nmax=10,
        tmax=tmax,
        mutationdist=:poisson,
        birthrate=1,
        deathrate=0.0,
        clonalmutations=0,
        μ=1
    )
    simulation = runsimulation(input, rng)
    @test length(allcells(simulation.output)) == 10
    @test age(simulation) <= tmax
end

@testset "updates" begin
    rng = MersenneTwister(12)
    cells = Vector{Union{BinaryNode{SimpleTreeCell}, Nothing}}([
        BinaryNode(SimpleTreeCell(id=1, birthtime=0.3, mutations=0, clonetype=1)),
        BinaryNode(SimpleTreeCell(id=2, birthtime=0.3, mutations=0, clonetype=1)),
        BinaryNode(SimpleTreeCell(id=3, birthtime=0.4, mutations=0, clonetype=1)),
    ])
    population = Population(cells, 1.0, 0.1, 1.0, 0.1)
    BirthDeathMutation.celldivision!(population, population.subclones, 1, 0.5, 5, [1], [:fixed], rng)
    @test length(allcells(population)) == 4
    @test getsubclonesizes(population) == [4]
    BirthDeathMutation.cellmutation!(population, population.subclones, 0.5, population.cells[1], 0.5)
    @test getclonetype(population.cells[1]) == 2
    @test length(population.subclones) == 2
    @test population.subclones[2].birthrate ≈ 1.5
    @test population.subclones[2].deathrate ≈ 0.1
    @test population.subclones[2].moranrate ≈ 1.5
    @test population.subclones[2].asymmetricrate ≈ 0.15
    @test getsubclonesizes(population) == [3, 1]
    n = length(allcells(population))
    BirthDeathMutation.moranupdate!(
        population,
        SelectionPredefined(Float64[0.5], Float64[0.5]),
        BirthDeathMutation.getmoranrates(population.subclones),
        maximum(BirthDeathMutation.getmoranrates(population.subclones)),
        n, 7, 2, 2, 0.51, [1], [:fixed], false, rng
    )
    @test length(allcells(population)) == n
end

@testset "rates" begin
    subclone = Subclone(1, 0, 0.0, 1, 2.0, 1.0, 1.0, 2.0)
    @test BirthDeathMutation.getwildtyperates([subclone]) == (birthrate=2.0, deathrate=1.0, moranrate=1.0, asymmetricrate=2.0)
    newrates = BirthDeathMutation.get_newsubclone_rates(BirthDeathMutation.getwildtyperates([subclone]), 0.5)
    @test newrates == (birthrate=3.0, deathrate=1.0, moranrate=1.5, asymmetricrate=3.0)
end

@testset "selection runsimulation" begin
    rng = MersenneTwister(100)
    input = BranchingInput(
        Nmax=10,
        mutationdist=:poisson,
        birthrate=1,
        deathrate=0.0,
        clonalmutations=0,
        μ=1
    )
    selection = SelectionPredefined([0.5], [3])
    simulation = runsimulation(input, selection, rng)
    @test length(allcells(simulation.output)) == 10
    @test sum(getsubclonesizes(simulation)) == 10
    @test getsubclonesizes(simulation) == counts(getclonetype.(allcells(simulation.output)), 1:length(simulation.output.subclones))

    tmax = 10
    input = MoranInput(
        N=10,
        tmax=tmax,
        mutationdist=:poisson,
        moranrate=1.0,
        clonalmutations=0,
        μ=1
    )
    simulation = runsimulation(input, rng)
    @test length(allcells(simulation.output)) == 10
    @test age(simulation) == simulation.output.t
    @test age(simulation) <= tmax

    tmax = 10
    input = BranchingMoranInput(
        Nmax=10,
        tmax=tmax,
        mutationdist=:poisson,
        birthrate=10,
        deathrate=0.0,
        clonalmutations=0,
        μ=1
    )
    selection = SelectionPredefined([0.5, 1.0], [0.1, 7.0])
    simulation = runsimulation(input, selection, rng)
    @test length(allcells(simulation.output)) == 10
    @test sum(getsubclonesizes(simulation)) == 10
    subclone_by_cell = getclonetype.(allcells(simulation.output))
    @test getsubclonesizes(simulation) == counts(subclone_by_cell, 1:length(simulation.output.subclones))
    @test age(simulation) <= tmax
    @test simulation.output.subclones[2].mutationtime ≈ 0.1 atol=0.1
    @test simulation.output.subclones[3].mutationtime ≈ 7.0 atol=1
end
