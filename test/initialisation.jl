@testset "population initialization" begin
    # Single cell, no initial mutations
    pop = initialize_population(1)
    @test length(allcells(pop)) == 1
    @test pop.t == 0.0
    @test allcells(pop)[1].data.mutations == 0
    @test allcells(pop)[1].data.clonetype == 1
    @test length(pop.subclones) == 1
    @test pop.subclones[1].s == 0.0

    # N cells with initial mutations
    pop = initialize_population(10; n_mutations=5)
    @test length(allcells(pop)) == 10
    @test all(cell.data.mutations == 5 for cell in allcells(pop))
    @test all(cell.data.clonetype == 1 for cell in allcells(pop))
    @test pop.subclones[1].size == 10

    # Clone-spec: two clones
    pop = initialize_population([
        (n=8, s=0.0, mutations=10),
        (n=2, s=0.2, mutations=15),
    ])
    @test length(allcells(pop)) == 10
    @test length(pop.subclones) == 2
    @test pop.subclones[1].s == 0.0
    @test pop.subclones[2].s == 0.2
    @test pop.subclones[1].size == 8
    @test pop.subclones[2].size == 2
    @test pop.t == 0.0

    cells = allcells(pop)
    @test all(cells[i].data.clonetype == 1 for i in 1:8)
    @test all(cells[i].data.clonetype == 2 for i in 9:10)
    @test all(cells[i].data.mutations == 10 for i in 1:8)
    @test all(cells[i].data.mutations == 15 for i in 9:10)
end
