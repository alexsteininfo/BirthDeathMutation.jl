@testset "tree mutations" begin
    rng = MersenneTwister(12)

    # Build a population via celldivision! to get a tree with known structure,
    # then verify final_timedep_mutations! applies time-dependent mutations correctly.
    input = BranchingInput(
        clonalmutations=0,
        Nmax=1,
        birthrate=1,
        deathrate=0,
        μ=[1, 3],
        mutationdist=[:fixed, :fixedtimedep]
    )
    population = initialize_population(input)
    subclones = population.subclones
    population, subclones, _ = BirthDeathMutation.celldivision!(
        population, subclones, 1, 1.0, 2, input.μ, input.mutationdist, rng
    )
    population, subclones, _ = BirthDeathMutation.celldivision!(
        population, subclones, 1, 3.0, 4, input.μ, input.mutationdist, rng
    )
    # cells[1]: born at t=1, last updated at t=3 (after 2nd division)
    # cells[2]: born at t=1, last updated at t=1
    # cells[3]: born at t=3, last updated at t=3
    BirthDeathMutation.final_timedep_mutations!(
        population, input.μ, input.mutationdist, rng; tend=4.0
    )
    # cells[1]: Δt = 4.0 - 3.0 = 1.0, adds round(3*1.0) = 3 → was 1, now 4
    # cells[2]: Δt = 4.0 - 1.0 = 3.0, adds round(3*3.0) = 9 → was 1, now 10
    # cells[3]: Δt = 4.0 - 3.0 = 1.0, adds round(3*1.0) = 3 → was 1, now 4
    @test population.cells[1].data.mutations == 4
    @test population.cells[2].data.mutations == 10
    @test population.cells[3].data.mutations == 4

    # changemutations! reassigns mutation counts to all nodes in the tree
    root = getsingleroot(allcells(population))
    BirthDeathMutation.changemutations!(root, 2, :fixed, 4.0, rng)
    for cellnode in AbstractTrees.PreOrderDFS(root)
        @test cellnode.data.mutations == 2
    end
end
