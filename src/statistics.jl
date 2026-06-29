"""
    mutations_per_cell(population::Population) -> Vector{Int64}

Calculate the total number of mutations per alive cell (summing along each cell's
lineage back to the root).
"""
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

"""
    average_mutations(population::Population) -> Float64

Calculate the mean number of mutations averaged over all alive cells.
"""
function average_mutations(population::Population)
    return mean(mutations_per_cell(population))
end

"""
    clonal_mutations(population::Population) -> Int64

Return the number of clonal mutations (mutations present in every alive cell),
computed as the total mutations at the MRCA node.
"""
function clonal_mutations(population::Population)
    MRCA = findMRCA(population)
    return all_cell_mutations(MRCA)
end

"""
    all_cell_mutations(node::BinaryNode) -> Int64

Return the total number of mutations for the cell at `node`, including inherited ones.
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
    pairwise_differences(population::Population[, idx]) -> Dict{Int64,Int64}
    pairwise_differences(cells::Vector{<:BinaryNode}[, idx]) -> Dict{Int64,Int64}

Calculate the number of pairwise differences between every pair of alive cells and
return as a dictionary (number of differences => frequency).
"""
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
    pairwisedistance(node1::BinaryNode, node2::BinaryNode) -> Int64

Return the number of mutational differences between two cells.
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
    pairwisedistances(population::Population[, idx]) -> Vector{Int64}

Return a vector of pairwise distances between every pair of alive cells.
"""
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

# Internal: used by coalescence_times
function _time_to_MRCA(node1, node2, t)
    if node1.data.birthtime > node2.data.birthtime
        node1, node2 = node2, node1
    end
    isdefined(node1, :parent) || return t - endtime(node1)
    isdefined(node2, :parent) || return t - endtime(node2)
    if node1.parent == node2.parent
        return t - node1.data.birthtime
    else
        return _time_to_MRCA(node1, node2.parent, t)
    end
end

"""
    coalescence_times(root[, idx]; t=nothing) -> Vector{Float64}
    coalescence_times(population::Population[, idx]; t=nothing) -> Vector{Float64}

Compute the coalescence time (time to MRCA) for every pair of alive cells.
"""
function coalescence_times(root::BinaryNode, idx=nothing; t=nothing)
    t = isnothing(t) ? age(root) : t
    coaltimes = Float64[]
    alivecells = getalivecells(root)
    alivecells = isnothing(idx) ? alivecells : alivecells[idx]
    while length(alivecells) > 1
        cellnode1 = popfirst!(alivecells)
        for cellnode2 in alivecells
            push!(coaltimes, _time_to_MRCA(cellnode1, cellnode2, t))
        end
    end
    return coaltimes
end

function coalescence_times(population::Population, idx=nothing; t=nothing)
    t = isnothing(t) ? age(population) : t
    coaltimes = Float64[]
    alivecells = allcells(population)
    alivecells = isnothing(idx) ? alivecells : alivecells[idx]
    while length(alivecells) > 1
        cellnode1 = popfirst!(alivecells)
        for cellnode2 in alivecells
            push!(coaltimes, _time_to_MRCA(cellnode1, cellnode2, t))
        end
    end
    return coaltimes
end

getsubclonesizes(subclones::Vector{Subclone}) = map(x -> length(x), subclones)
getsubclonesizes(population::AbstractPopulation) = getsubclonesizes(population.subclones)

"""
    sitefrequencyspectrum(population::Population) -> Vector{Int64}

Compute the site-frequency spectrum from the neutral mutation tree.
`sfs[k]` = number of neutral mutations present in exactly `k` alive cells.
"""
function sitefrequencyspectrum(population::Population)
    sfs  = zeros(Int64, popsize(population))
    root = getsingleroot(allcells(population))
    _sfs_fill!(root, sfs)
    return sfs
end

function _sfs_fill!(node::BinaryNode, sfs::Vector{Int64})
    if isnothing(node.left) && isnothing(node.right)
        sfs[1] += node.data.mutations
        return 1
    end
    count = 0
    isnothing(node.left)  || (count += _sfs_fill!(node.left,  sfs))
    isnothing(node.right) || (count += _sfs_fill!(node.right, sfs))
    count > 0 && (sfs[count] += node.data.mutations)
    return count
end
