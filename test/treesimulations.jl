function make_tree()
    root = BinaryNode(SimpleTreeCell(id=1, birthtime=0.0, mutations=0))
    leftchild!(root, SimpleTreeCell(id=2, birthtime=1.5355542835848743, mutations=5))
    rightchild!(root, SimpleTreeCell(id=3, birthtime=1.5355542835848743, mutations=10))
    leftchild!(root.left, SimpleTreeCell(id=4, birthtime=1.6919799338516708, mutations=13))
    rightchild!(root.left, SimpleTreeCell(id=5, birthtime=1.6919799338516708, mutations=17))
    return root
end

function make_tree2()
    root = BinaryNode(SimpleTreeCell(id=1, birthtime=2.0, mutations=0))
    leftchild!(root, SimpleTreeCell(id=1, birthtime=3.0, mutations=5))
    rightchild!(root, SimpleTreeCell(id=1, birthtime=3.0, mutations=11))
    return root
end

@testset "tree branching" begin
    rng = MersenneTwister(12)

    @testset "division" begin
        input = BranchingInput(clonalmutations=0, Nmax=10, birthrate=1, deathrate=0)
        population = initialize_population(input)
        subclones = population.subclones
        root = getsingleroot(allcells(population))
        population, subclones, nextID = BirthDeathMutation.celldivision!(
            population, subclones, 1, 1.1, 2, [10], [:fixedtimedep], rng
        )
        @test nextID == 4
        @test length(allcells(population)) == 2
        @test !(root in allcells(population))
        @test root.data.mutations == 11
        for (i, cellnode) in enumerate(allcells(population))
            @test cellnode.parent == root
            @test cellnode.data.birthtime ≈ 1.1 atol=0.001
            @test cellnode.data.mutations == 0
            @test cellnode.data.clonetype == root.data.clonetype
            @test cellnode.data.id == i + 1
        end
        BirthDeathMutation.celldivision!(
            population, subclones, 2, 1.1, 3, [10], [:fixedtimedep], rng; nchildcells=1
        )
        @test length(allcells(population)) == 2
        @test !(root in allcells(population))
        @test population.cells[end] == root.right.left
        @test isnothing(root.right.right)
        @test population.cells[end].data.id == 3
    end

    @testset "death simple" begin
        input = BranchingInput(clonalmutations=0, Nmax=1, birthrate=1, deathrate=0)
        population = initialize_population(input)
        subclones = population.subclones
        population, subclones, _ = BirthDeathMutation.celldivision!(
            population, subclones, 1, 1.1, 2, [10], [:fixedtimedep], rng
        )
        @test length(subclones[1]) == 2
        BirthDeathMutation.cellmutation!(population, subclones, 0.5, population.cells[1], 1.1)
        @test length(subclones) == 2
        @test length(subclones[1]) == 1
        @test length(subclones[2]) == 1
        @test subclones[2].subcloneid == 2
        @test subclones[2].parentid == 1
        @test subclones[2].mutationtime == 1.1
        @test subclones[2].size == 1
        @test subclones[2].birthrate == 1.5
        @test subclones[2].deathrate == 0.0
        @test subclones[2].moranrate == 0.0
        @test subclones[2].asymmetricrate == 0.0

        BirthDeathMutation.celldeath!(population, subclones, 1, 1.5, [10], [:fixedtimedep], rng)
        @test length(allcells(population)) == 1
        BirthDeathMutation.celldivision!(
            population, subclones, 1, 1.1, 2, [10], [:fixedtimedep], rng
        )
        idx = 2
        deadcellnode = population.cells[idx]
        BirthDeathMutation.celldeath!(population, subclones, idx, 2.0, [10], [:fixedtimedep], rng)
        @test !(deadcellnode in allcells(population))
    end

    @testset "mutations" begin
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
        @test population.cells[1].data.mutations == 0
        population, subclones, _ = BirthDeathMutation.celldivision!(
            population, subclones, 1, 1.0, 2, input.μ, input.mutationdist, rng
        )
        @test population.cells[1].data.mutations == 1
        @test population.cells[2].data.mutations == 1
        @test population.cells[1].parent.data.mutations == 3
        population, subclones, _ = BirthDeathMutation.celldivision!(
            population, subclones, 1, 3.0, 4, input.μ, input.mutationdist, rng
        )
        @test population.cells[1].data.mutations == 1
        @test population.cells[2].data.mutations == 1
        @test population.cells[3].data.mutations == 1
        @test population.cells[1].parent.data.mutations == 7
        BirthDeathMutation.final_timedep_mutations!(
            population, input.μ, input.mutationdist, rng; tend=4.0
        )
        @test population.cells[1].data.mutations == 4
        @test population.cells[2].data.mutations == 10
        @test population.cells[3].data.mutations == 4
    end
end

tree_root = make_tree()

@testset "basic tree" begin
    root = tree_root
    @test AbstractTrees.nextsibling(root.left.left) == root.left.right
    @test isnothing(AbstractTrees.nextsibling(root.left.right))
    @test AbstractTrees.prevsibling(root.left.right) == root.left.left
    @test isnothing(AbstractTrees.prevsibling(root.left.left))
    @test isnothing(AbstractTrees.nextsibling(root))
end

@testset "prune tree" begin
    root = make_tree()
    BirthDeathMutation.prune_tree!(root.left.left)
    @test !intree(root.left.left, root)
    BirthDeathMutation.prune_tree!(root.left.right)
    @test !intree(root.left, root)
end

@testset "tree statistics" begin
    root = tree_root
    @test mutations_per_cell(root, includeclonal=true) == [18, 22, 10]
    @test celllifetimes(root, excludeliving=true) ≈ [1.53555, 0.15643] atol=0.01
    @test celllifetimes(root, excludeliving=false) ≈ [1.53555, 0.15643, 0.0, 0.0, 0.15643] atol=0.01
    @test time_to_MRCA(root.left.left, root.right, 2.0) ≈ 2.0 - 1.5355542835848743
    @test time_to_MRCA(root.left.left, root.left.right, 2.0) ≈ 2.0 - 1.6919799338516708
    @test coalescence_times(root) ≈ [0.0, 1.6919799338516708 - 1.5355542835848743, 1.6919799338516708 - 1.5355542835848743]
    @test pairwisedistance(root.left.left, root.right) == 28
    @test pairwisedistance(root.left.left, root.left.right) == 30
    @test pairwise_differences(getalivecells(root)) == Dict(28=>1, 30=>1, 32=>1)
    @test Set(getalivecells(root)) == Set([root.right, root.left.left, root.left.right])
    @test endtime(root.left.left) === nothing
    @test endtime(root) ≈ 1.5355542835848743 atol=1e-6
end

tree_root2 = make_tree2()
tree_roots = vcat(tree_root, tree_root2)
alivecells_2roots = [tree_root.right, tree_root.left.left, tree_root.left.right, tree_root2.left, tree_root2.right]

@testset "multiple roots" begin
    @test getroot([tree_root2.left, tree_root2.right]) == [tree_root2]
    @test Set(getroot(alivecells_2roots)) == Set(tree_roots)
    @test Set(getalivecells(tree_roots)) == Set(alivecells_2roots)
    @test isnothing(getsingleroot(alivecells_2roots))
    @test getsingleroot([tree_root2.left, tree_root2.right]) == tree_root2
end

@testset "MRCA" begin
    root = make_tree()
    alivecells = collect(Leaves(root))
    @test findMRCA(alivecells) == root
    @test findMRCA(root.left.left, root.left.right) == root.left
    @test findMRCA(root.left.left, root.right) == root
end

@testset "updates" begin
    rng = MersenneTwister(12)
    t2 = make_tree2()
    cells = Vector{Union{BinaryNode{SimpleTreeCell}, Nothing}}(
        [cellnode for cellnode in Leaves(t2)]
    )
    population = Population(deepcopy(cells), 1.0, 0.1, 1.0, 0.1)
    BirthDeathMutation.celldivision!(population, population.subclones, 1, 3.5, 3, [1], [:fixed], rng)
    @test length(allcells(population)) == 3
    @test getsubclonesizes(population) == [3]
    BirthDeathMutation.cellmutation!(population, population.subclones, 0.5, population.cells[1], 3.5)
    @test getclonetype(population.cells[1]) == 2
    @test length(population.subclones) == 2
    @test population.subclones[2].birthrate ≈ 1.5
    @test population.subclones[2].deathrate ≈ 0.1
    @test population.subclones[2].moranrate ≈ 1.5
    @test population.subclones[2].asymmetricrate ≈ 0.15
    @test getsubclonesizes(population) == [2, 1]
    n = length(allcells(population))
    BirthDeathMutation.moranupdate!(
        population,
        SelectionPredefined([3.5], [0.5]),
        BirthDeathMutation.getmoranrates(population.subclones),
        maximum(BirthDeathMutation.getmoranrates(population.subclones)),
        n, 7, 2, 2, 4.0, [1], [:fixed], false, rng
    )
    @test length(allcells(population)) == n
end

@testset "selection tree runsimulation" begin
    rng = MersenneTwister(100)
    input = BranchingInput(
        Nmax=10,
        mutationdist=:poisson,
        birthrate=1,
        deathrate=0.0,
        clonalmutations=0,
        μ=1,
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
    selection = SelectionPredefined([0.5, 1.0], [0.1, 8])
    simulation = runsimulation(input, selection, rng)
    @test length(allcells(simulation.output)) == 10
    @test sum(getsubclonesizes(simulation)) == 10
    subclone_by_cell = getclonetype.(allcells(simulation.output))
    @test getsubclonesizes(simulation) == counts(subclone_by_cell, 1:length(simulation.output.subclones))
    @test age(simulation) <= tmax
    @test simulation.output.subclones[2].mutationtime ≈ 0.1 atol=0.1
    @test simulation.output.subclones[3].mutationtime ≈ 8 atol=1
end
