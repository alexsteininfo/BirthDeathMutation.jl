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
        pop = initialize_population(1)
        subclones = pop.subclones
        root = getsingleroot(allcells(pop))
        pop, subclones, nextID = BirthDeathMutation.celldivision!(
            pop, subclones, 1, 1.1, 2, 1.0, :fixed, rng
        )
        @test nextID == 4
        @test length(allcells(pop)) == 2
        @test !(root in allcells(pop))
        @test root.data.mutations == 0  # parent gets no mutations with :fixed μ=1 in children
        for (i, cellnode) in enumerate(allcells(pop))
            @test cellnode.parent == root
            @test cellnode.data.birthtime ≈ 1.1 atol=0.001
            @test cellnode.data.mutations == 1  # :fixed, μ=1
            @test cellnode.data.clonetype == root.data.clonetype
            @test cellnode.data.id == i + 1
        end
        BirthDeathMutation.celldivision!(
            pop, subclones, 2, 1.1, 3, 0.0, :fixed, rng; nchildcells=1
        )
        @test length(allcells(pop)) == 2
        @test !(root in allcells(pop))
        @test pop.cells[end] == root.right.left
        @test isnothing(root.right.right)
        @test pop.cells[end].data.id == 3
    end

    @testset "death" begin
        pop = initialize_population(1)
        subclones = pop.subclones
        pop, subclones, _ = BirthDeathMutation.celldivision!(
            pop, subclones, 1, 1.1, 2, 0.0, :fixed, rng
        )
        @test length(subclones[1]) == 2
        BirthDeathMutation.cellmutation!(pop, subclones, 0.5, pop.cells[1], 1.1)
        @test length(subclones) == 2
        @test length(subclones[1]) == 1
        @test length(subclones[2]) == 1
        @test subclones[2].subcloneid == 2
        @test subclones[2].parentid == 1
        @test subclones[2].mutationtime == 1.1
        @test subclones[2].size == 1
        @test subclones[2].s == 0.5

        BirthDeathMutation.celldeath!(pop, subclones, 1, 1.5)
        @test length(allcells(pop)) == 1
        BirthDeathMutation.celldivision!(
            pop, subclones, 1, 1.1, 2, 0.0, :fixed, rng
        )
        idx = 2
        deadcellnode = pop.cells[idx]
        BirthDeathMutation.celldeath!(pop, subclones, idx, 2.0)
        @test !(deadcellnode in allcells(pop))
    end

    @testset "mutations" begin
        pop = initialize_population(1)
        subclones = pop.subclones
        @test pop.cells[1].data.mutations == 0
        pop, subclones, _ = BirthDeathMutation.celldivision!(
            pop, subclones, 1, 1.0, 2, 3.0, :fixed, rng
        )
        @test pop.cells[1].data.mutations == 3
        @test pop.cells[2].data.mutations == 3
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

@testset "selection tree simulate!" begin
    rng = MersenneTwister(100)
    Nmax = 10
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

    # Two selection events in a BirthDeathBlock
    rng = MersenneTwister(100)
    Nmax = 20
    sel = SelectionPredefined([0.5, 1.0], [0.1, 2.0])
    block2 = BirthDeathBlock(
        birthrate    = (s, N) -> 10.0 * (1 + s),
        deathrate    = (s, N) -> 0.0,
        stopfunction = pop -> popsize(pop) >= Nmax,
        selection    = sel,
        μ = 1.0,
        mutationdist = :poisson,
    )
    pop = initialize_population(1)
    pop = simulate!(pop, block2, rng)
    @test length(allcells(pop)) == Nmax
    @test sum(getsubclonesizes(pop)) == Nmax
    # mutant must appear at or after the specified time
    @test pop.subclones[2].mutationtime >= sel.mutant_time[1]
end

@testset "changemutations!" begin
    rng = MersenneTwister(12)
    root = make_tree()
    BirthDeathMutation.changemutations!(root, 2, :fixed, rng)
    for cellnode in AbstractTrees.PreOrderDFS(root)
        @test cellnode.data.mutations == 2
    end
end
