
"""
    killcell!(alivecells, deadcellidx, args...)

Kill `SimpleTreeCell` at index `deadcellidx` in `alivecells` by pruning all references
to it from the tree.
"""
function killcell!(alivecells, deadcellidx::Integer, args...)
    prune_tree!(alivecells[deadcellidx])
    return alivecells
end

"""
    killcells!(alivecells, deadcellvector, args...)

Kill multiple tree cells.
"""
function killcells!(alivecells, deadcellvector, args...)
    for deadcellidx in deadcellvector
        killcell!(alivecells, deadcellidx, args...)
    end
    return alivecells
end

function killallcells!(alivecells, args...)
    killcells!(alivecells, 1:lastindex(alivecells), args...)
    return alivecells
end

"""
    prune_tree!(cellnode)

Remove `cellnode` from tree, and remove any ancestor node that has no remaining children.
"""
function prune_tree!(cellnode)
    while true
        parent = cellnode.parent
        if isnothing(parent)
            return
        else
            cellnode.parent = nothing
            if parent.left == cellnode
                parent.left = nothing
            elseif parent.right == cellnode
                parent.right = nothing
            else
                error("dead cell is neither left nor right child of parent")
            end
            if isnothing(parent.left) && isnothing(parent.right)
                cellnode = parent
            else
                return
            end
        end
    end
end

"""
    changemutations!(root::BinaryNode, μ, mutationdist, rng, clonalmutations=0)

Reassign mutations to every node in the phylogeny rooted at `root`.
"""
function changemutations!(root::BinaryNode, μ, mutationdist, rng, clonalmutations=0)
    for cellnode in PreOrderDFS(root)
        cellnode.data.mutations = numbernewmutations(rng, mutationdist, Float64(μ))
    end
    root.data.mutations += clonalmutations
end

"""
    endtime(cellnode::BinaryNode)

Return the time at which the cell divided. If the cell is still alive, return `nothing`.
"""
function endtime(cellnode::BinaryNode)
    if haschildren(cellnode)
        return cellnode.left.data.birthtime
    else
        return nothing
    end
end

"""
    celllifetime(cellnode::BinaryNode, [tmax])

Compute the lifetime of a cell. If it hasn't divided yet, use `tmax` as the end time.
"""
function celllifetime(cellnode::BinaryNode, tmax=nothing)
    if haschildren(cellnode)
        return cellnode.left.data.birthtime - cellnode.data.birthtime
    else
        tmax = isnothing(tmax) ? age(getroot(cellnode)) : tmax
        return tmax - cellnode.data.birthtime
    end
end

"""
    celllifetimes(root; excludeliving=true)

Compute the lifetime of each cell in the phylogeny, excluding currently alive cells by
default.
"""
function celllifetimes(root; excludeliving=true)
    lifetimes = Float64[]
    if excludeliving
        for cellnode in PreOrderDFS(root)
            if haschildren(cellnode)
                push!(lifetimes, cellnode.left.data.birthtime - cellnode.data.birthtime)
            end
        end
    else
        popage = age(root)
        for cellnode in PreOrderDFS(root)
            push!(lifetimes, celllifetime(cellnode, popage))
        end
    end
    return lifetimes
end

"""
    age(root::BinaryNode)

Compute the age of the population as the birthtime of the most recently born leaf cell.
"""
function age(root::BinaryNode)
    age = 0
    for cellnode in Leaves(root)
        if cellnode.data.birthtime > age
            age = cellnode.data.birthtime
        end
    end
    return age
end

getalivecells(root::BinaryNode) =
    [cellnode for cellnode in Leaves(root) if isalive(cellnode.data)]

getalivecells(roots::Vector{BinaryNode{T}}) where T =
    [cellnode for root in roots for cellnode in Leaves(root) if isalive(cellnode.data)]

isalive(cellnode::BinaryNode{T}) where T = isalive(cellnode.data)
isalive(cell::SimpleTreeCell) = true
isalive(::Nothing) = false

id(cellnode::BinaryNode{<:AbstractTreeCell}) = cellnode.data.id

"""
    asroot!(node)

Transform `node` into a root by setting its `parent` field to `nothing`. Returns `node`
and the original parent.
"""
function asroot!(node)
    parent = node.parent
    node.parent = nothing
    return node, parent
end

"""
    cell_subset_size(node, cells)

Count the alive leaf nodes under `node` that are present in `cells`.
"""
function cell_subset_size(node, cells)
    if node in cells
        if isalive(node) return 1 else 0 end
    end
    root, parent = asroot!(node)
    count = 0
    for leaf in Leaves(root)
        if leaf in cells && isalive(leaf)
            count += 1
        end
    end
    root.parent = parent
    return count
end
