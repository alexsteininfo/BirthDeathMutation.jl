"""
    initialize_population(input; rng=Random.GLOBAL_RNG)

Create the initial `Population` for the given `input`.
"""
function initialize_population(input::SinglelevelInput)
    N = getNinit(input)
    cells = Vector{Union{BinaryNode{SimpleTreeCell}, Nothing}}(create_cells(input.clonalmutations, N))
    birthrate, deathrate, moranrate, asymmetricrate = getinputrates(input)
    return Population(cells, birthrate, deathrate, moranrate, asymmetricrate)
end

getNinit(input::Union{BranchingInput, BranchingMoranInput}) = 1
getNinit(input::MoranInput) = input.N

function newcell(id, mutations)
    return BinaryNode{SimpleTreeCell}(SimpleTreeCell(;id, mutations))
end

function create_cells(initialmutations, N=1)
    return [newcell(id, initialmutations) for id in 1:N]
end
