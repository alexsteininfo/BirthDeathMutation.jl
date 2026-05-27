@testset "population initialization" begin
    # BranchingInput: starts with 1 cell, clonal mutations assigned to root
    input = BranchingInput(Nmax=100, birthrate=1.0, deathrate=0.0, clonalmutations=50, μ=1)
    population = initialize_population(input)
    @test length(allcells(population)) == 1
    @test population.t == 0.0
    @test allcells(population)[1].data.mutations == 50
    @test allcells(population)[1].data.clonetype == 1
    @test length(population.subclones) == 1
    @test BirthDeathMutation.getwildtyperates(population) == (
        birthrate=1.0, deathrate=0.0, moranrate=0.0, asymmetricrate=0.0
    )

    # MoranInput: starts with N identical cells
    input = MoranInput(N=10, tmax=10.0, moranrate=2.0, clonalmutations=5, μ=1)
    population = initialize_population(input)
    @test length(allcells(population)) == 10
    @test population.t == 0.0
    @test all(cell.data.mutations == 5 for cell in allcells(population))
    @test all(cell.data.clonetype == 1 for cell in allcells(population))
    @test BirthDeathMutation.getwildtyperates(population) == (
        birthrate=0.0, deathrate=0.0, moranrate=2.0, asymmetricrate=0.0
    )

    # BranchingMoranInput: starts with 1 cell
    input = BranchingMoranInput(
        Nmax=50, tmax=10.0, birthrate=2.0, deathrate=0.5, clonalmutations=0, μ=1
    )
    population = initialize_population(input)
    @test length(allcells(population)) == 1
    @test population.t == 0.0
    @test BirthDeathMutation.getwildtyperates(population) == (
        birthrate=2.0, deathrate=0.5, moranrate=1.0, asymmetricrate=0.0
    )
end
