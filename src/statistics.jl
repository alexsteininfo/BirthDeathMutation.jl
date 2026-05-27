
"""
    mutations_per_cell(simulation::Simulation)
    mutations_per_cell(population::Population)

Calculate the number of mutations per alive cell (summing mutations along each cell's
lineage back to the root).
"""
mutations_per_cell(simulation::Simulation) = mutations_per_cell(simulation.output)

function mutations_per_cell(root::BinaryNode{T}; includeclonal=false) where T <: AbstractTreeCell
    mutspercell = Int64[]
    for cellnode in Leaves(root)
        if isalive(cellnode.data)
            mutations = cellnode.data.mutations
            node = cellnode
            while true
                if !AbstractTrees.isroot(node) && (node != root || includeclonal)
                    node = node.parent
                    mutations += node.data.mutations
                else
                    break
                end
            end
            push!(mutspercell, mutations)
        end
    end
    return mutspercell
end

function mutations_per_cell(population::Population)
    mutspercell = Int64[]
    for cellnode in allcells(population)
        mutations = cellnode.data.mutations
        node = cellnode
        while !isnothing(node.parent)
            node = node.parent
            mutations += node.data.mutations
        end
        push!(mutspercell, mutations)
    end
    return mutspercell
end

"""
    average_mutations(simulation::Simulation)
    average_mutations(population::Population)

Calculate the mean number of mutations averaged over all alive cells.
"""
average_mutations(simulation::Simulation) = average_mutations(simulation.output)

function average_mutations(population::Population)
    return mean(mutations_per_cell(population))
end

"""
    clonal_mutations(simulation::Simulation)
    clonal_mutations(population::Population)

Return the number of clonal mutations (mutations present in every alive cell), computed
as the total mutations at the MRCA node.
"""
function clonal_mutations(simulation::Simulation)
    return clonal_mutations(simulation.output)
end

function clonal_mutations(population::Population)
    MRCA = findMRCA(population)
    return all_cell_mutations(MRCA)
end

"""
    all_cell_mutations(node::BinaryNode)

Return the total number of mutations for the cell at `node`, including those inherited.
"""
function all_cell_mutations(node::BinaryNode)
    mutations = node.data.mutations
    while !isnothing(node.parent)
        node = node.parent
        mutations += node.data.mutations
    end
    return mutations
end

"""
    pairwise_differences(simulation::Simulation[, idx])
    pairwise_differences(population::Population[, idx])
    pairwise_differences(cells::Vector{<:BinaryNode}[, idx])

Calculate the number of pairwise differences between every pair of alive cells and return
as a dictionary (number differences => frequency). If `idx` is given, only include listed
cells.
"""
function pairwise_differences(simulation::Simulation, idx=nothing)
    return pairwise_differences(simulation.output, idx)
end

function pairwise_differences(population::Population, idx=nothing)
    return pairwise_differences(allcells(population), idx)
end

function pairwise_differences(cells::Vector{<:BinaryNode}, idx=nothing)
    if !isnothing(idx)
        cells = cells[idx]
    end
    n = length(cells)
    pfd_vec = Int64[]
    for i in 1:n
        for j in i+1:n
            push!(pfd_vec, pairwisedistance(cells[i], cells[j]))
        end
    end
    return countmap(pfd_vec)
end

"""
    pairwisedistance(node1::BinaryNode, node2::BinaryNode)

Return the number of pairwise differences between two cells.
"""
function pairwisedistance(cellnode1::BinaryNode, cellnode2::BinaryNode)
    cellnode1 == cellnode2 && return 0
    distance = 0
    while true
        if cellnode1.data.id > cellnode2.data.id
            cellnode1, cellnode2 = cellnode2, cellnode1
        end
        if (
            (cellnode1.parent == cellnode2.parent) ||
            (isnothing(cellnode1.parent) && isnothing(cellnode2.parent))
        )
            return distance + cellnode1.data.mutations + cellnode2.data.mutations
        elseif cellnode1 == cellnode2.parent
            return distance + cellnode2.data.mutations
        else
            distance += cellnode2.data.mutations
            cellnode2 = cellnode2.parent
        end
    end
end

"""
    pairwisedistances(simulation::Simulation[, idx])
    pairwisedistances(population::Population[, idx])

Return a vector of pairwise distances between every pair of alive cells. If `idx` is
given, only include listed cells.
"""
function pairwisedistances(simulation::Simulation, idx=nothing)
    return pairwisedistances(simulation.output, idx)
end

function pairwisedistances(population::Population, idx=nothing)
    cells = allcells(population)
    if !isnothing(idx)
        cells = cells[idx]
    end
    n = length(cells)
    pfd_vec = Int64[]
    for i in 1:n
        for j in i+1:n
            push!(pfd_vec, pairwisedistance(cells[i], cells[j]))
        end
    end
    return pfd_vec
end

"""
    pairwise_fixed_differences(simulation::Simulation)
    pairwise_fixed_differences(population::Population)

Alias for `pairwise_differences` — returns pairwise distances as a frequency dictionary.
"""
pairwise_fixed_differences(simulation::Simulation, idx=nothing) =
    pairwise_differences(simulation, idx)
pairwise_fixed_differences(population::Population, idx=nothing) =
    pairwise_differences(population, idx)

"""
    pairwise_fixed_differences_clonal(simulation::Simulation)
    pairwise_fixed_differences_clonal(population::Population)

Return pairwise distance frequency dict and clonal mutations frequency dict.
"""
function pairwise_fixed_differences_clonal(simulation::Simulation, idx=nothing)
    return pairwise_fixed_differences_clonal(simulation.output, idx)
end

function pairwise_fixed_differences_clonal(population::Population, idx=nothing)
    pfdvec = pairwisedistances(population, idx)
    clonalmuts = clonal_mutations(population)
    return countmap(pfdvec), countmap([clonalmuts])
end

"""
    pairwise_fixed_differences_matrix(simulation::Simulation; diagonals=false)
    pairwise_fixed_differences_matrix(population::Population; diagonals=false)

Return an n×n matrix of pairwise distances between all alive cells. If `diagonals=true`,
diagonal entries contain total mutations per cell.
"""
function pairwise_fixed_differences_matrix(simulation::Simulation, idx=nothing; diagonals=false)
    return pairwise_fixed_differences_matrix(simulation.output, idx; diagonals)
end

function pairwise_fixed_differences_matrix(population::Population, idx=nothing; diagonals=false)
    cells = allcells(population)
    if !isnothing(idx)
        cells = cells[idx]
    end
    n = length(cells)
    pfd = zeros(Int64, n, n)
    for i in 1:n
        if diagonals pfd[i,i] = all_cell_mutations(cells[i]) end
        for j in i+1:n
            pfd[j, i] = pairwisedistance(cells[i], cells[j])
        end
    end
    return pfd
end

"""
    pairwise_fixed_differences_statistics(simulation::Simulation[, samplesize, idx])
    pairwise_fixed_differences_statistics(population::Population[, idx])

Calculate mean and variance of pairwise distances and of clonal mutations.
"""
function pairwise_fixed_differences_statistics(simulation::Simulation, idx=nothing)
    pairwise_fixed_differences_statistics(simulation.output, idx)
end

function pairwise_fixed_differences_statistics(simulation::Simulation, samplesize::Integer, rng)
    pairwise_fixed_differences_statistics(simulation.output, samplesize, rng)
end

function pairwise_fixed_differences_statistics(population::Population, samplesize::Int64, rng)
    idx = sample(rng, 1:length(allcells(population)), samplesize, replace=false)
    return pairwise_fixed_differences_statistics(population, idx)
end

function pairwise_fixed_differences_statistics(population::Population, idx=nothing)
    pfd = pairwisedistances(population, idx)
    clonalmuts = clonal_mutations(population)
    return mean(pfd), var(pfd), clonalmuts, 0.0
end

"""
    time_to_MRCA(node1::BinaryNode, node2::BinaryNode, t)

Computes the time elapsed since the MRCA of two cells divided, relative to time `t`.
"""
function time_to_MRCA(node1, node2, t)
    if node1.data.birthtime > node2.data.birthtime
        node1, node2 = node2, node1
    end
    isdefined(node1, :parent) || return t - endtime(node1)
    isdefined(node2, :parent) || return t - endtime(node2)
    if node1.parent == node2.parent
        return t - node1.data.birthtime
    else
         return time_to_MRCA(node1, node2.parent, t)
    end
end

"""
    coalescence_times(root, [idx]; t=nothing)

Computes the coalescence time (time to MRCA) for every pair of alive cells under `root`.
"""
function coalescence_times(root, idx=nothing; t=nothing)
    t = isnothing(t) ? age(root) : t
    coaltimes = Float64[]
    alivecells = getalivecells(root)
    alivecells = isnothing(idx) ? alivecells : alivecells[idx]
    while length(alivecells) > 1
        cellnode1 = popfirst!(alivecells)
        for cellnode2 in alivecells
            push!(coaltimes, time_to_MRCA(cellnode1, cellnode2, t))
        end
    end
    return coaltimes
end

getsubclonesizes(subclones::Vector{Subclone}) = map(x -> length(x), subclones)
getsubclonesizes(population::AbstractPopulation) = getsubclonesizes(population.subclones)
getsubclonesizes(simulation::Simulation) = getsubclonesizes(simulation.output)
