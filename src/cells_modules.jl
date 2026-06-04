#region Define BinaryNode and methods for binary trees

"""
    BinaryNode{T}

Basic unit of a binary tree.

# Fields
- `data::T`
- `parent::Union{Nothing, BinaryNode{T}}`
- `left::Union{Nothing, BinaryNode{T}}`
- `right::Union{Nothing, BinaryNode{T}}`
"""
mutable struct BinaryNode{T}
    data::T
    parent::Union{Nothing, BinaryNode{T}}
    left::Union{Nothing, BinaryNode{T}}
    right::Union{Nothing, BinaryNode{T}}

    function BinaryNode{T}(data, parent=nothing, l=nothing, r=nothing) where T
        new{T}(data, parent, l, r)
    end
end
BinaryNode(data) = BinaryNode{typeof(data)}(data)

"""
    leftchild!(parent::BinaryNode, data)

Create a new `BinaryNode` from `data` and assign it to `parent.left`.

See also [`rightchild!`](@ref).
"""
function leftchild!(parent::BinaryNode, data)
    isnothing(parent.left) || error("left child is already assigned")
    node = typeof(parent)(data, parent)
    parent.left = node
end

"""
    rightchild!(parent::BinaryNode, data)

Create a new `BinaryNode` from `data` and assign it to `parent.right`.

See also [`leftchild!`](@ref).
"""
function rightchild!(parent::BinaryNode, data)
    isnothing(parent.right) || error("right child is already assigned")
    node = typeof(parent)(data, parent)
    parent.right = node
end

function AbstractTrees.children(node::BinaryNode)
    if isnothing(node.left) && isnothing(node.right)
        ()
    elseif isnothing(node.left) && !isnothing(node.right)
        (node.right,)
    elseif !isnothing(node.left) && isnothing(node.right)
        (node.left,)
    else
        (node.left, node.right)
    end
end

function AbstractTrees.nextsibling(child::BinaryNode)
    isnothing(child.parent) && return nothing
    p = child.parent
    if !isnothing(p.right)
        child === p.right && return nothing
        return p.right
    end
    return nothing
end

function AbstractTrees.prevsibling(child::BinaryNode)
    isnothing(child.parent) && return nothing
    p = child.parent
    if !isnothing(p.left)
        child === p.left && return nothing
        return p.left
    end
    return nothing
end

AbstractTrees.nodevalue(n::BinaryNode) = n.data
AbstractTrees.ParentLinks(::Type{<:BinaryNode}) = StoredParents()
AbstractTrees.parent(n::BinaryNode) = n.parent
AbstractTrees.NodeType(::Type{<:BinaryNode{T}}) where {T} = HasNodeType()
AbstractTrees.nodetype(::Type{<:BinaryNode{T}}) where {T} = BinaryNode{T}

Base.eltype(::Type{<:TreeIterator{BinaryNode{T}}}) where T = BinaryNode{T}
Base.IteratorEltype(::Type{<:TreeIterator{BinaryNode{T}}}) where T = Base.HasEltype()

AbstractTrees.printnode(io::IO, node::BinaryNode) = print(io, node.data)

"""
    popsize(root::BinaryNode)

Get the number of alive cells that are descendents of `root`.
"""
function popsize(root::BinaryNode)
    N = 0
    for l in Leaves(root)
        if isalive(nodevalue(l))
            N += 1
        end
    end
    return N
end

haschildren(node::BinaryNode) = length(children(node)) != 0

function Base.show(io::IO, node::BinaryNode)
    show(io, node.data)
end

function Base.show(io::IO, nodevec::Vector{BinaryNode{T}}) where T
    println(io, "$(length(nodevec))-element Vector{BinaryNode{$T}}:")
    for node in nodevec
        show(io, node)
        print(io, "\n")
    end
end
#endregion

#region Define cell types
abstract type AbstractCell end
abstract type AbstractTreeCell <: AbstractCell end

"""
    SimpleTreeCell

Represents a single cell that can be the `data` field of a `BinaryNode{SimpleTreeCell}`.
Dead cells are pruned from the tree structure.
"""
mutable struct SimpleTreeCell <: AbstractTreeCell
    id::Int64
    birthtime::Float64
    mutations::Int64
    clonetype::Int64
end

function SimpleTreeCell(; id=1, birthtime=0.0, mutations=0, clonetype=1)
    return SimpleTreeCell(id, birthtime, mutations, clonetype)
end

#endregion

"""
    Subclone

Defines subclone properties. `s` is the selection coefficient used in the block's
rate functions `b(s, N)` and `d(s, N)` (or `r(s, N)` for Moran blocks).
"""
@kwdef mutable struct Subclone
    subcloneid::Int64 = 1
    parentid::Int64 = 0
    mutationtime::Float64 = 0.0
    size::Int64 = 1
    s::Float64 = 0.0
    asymmetricrate::Float64 = 0.0
end

Base.length(subclone::Subclone) = subclone.size

getclonetype(cellnode::BinaryNode) = cellnode.data.clonetype
setclonetype!(cellnode::BinaryNode, newclonetype) = (cellnode.data.clonetype = newclonetype)

#region Methods for finding most recent common ancestor and tree roots

"""
    getroot(nodevec)

Get a list of unique roots for each tree that contains a node in `nodevec`.
This requires node to have the trait StoredParents.
"""
AbstractTrees.getroot(nodevec::Vector{BinaryNode{T}}) where T =
    collect(Set(getroot.(nodevec)))

"""
    getsingleroot(nodevec)

Find the unique root of a list of nodes `nodevec`. If there is no unique root return
`nothing`.
"""
function getsingleroot(nodevec::Vector{BinaryNode{T}}) where T
    roots = AbstractTrees.getroot(nodevec)
    if length(roots) == 1
        return roots[1]
    else
        return nothing
    end
end

"""
    findMRCA(population)
    findMRCA(nodes)
    findMRCA(node1, node2)

Find the most recent common ancestor of all nodes in `population.cells`, a list of nodes
`nodes`, or two nodes `node1` and `node2`.
"""
function findMRCA end

function findMRCA(node1, node2)
    node1 == node2 && return node1
    while true
        if  node1.data.id > node2.data.id
            node1, node2 = node2, node1
        end
        if node1.parent == node2.parent
            return node1.parent
        elseif isnothing(node1.parent) && isnothing(node2.parent)
            return nothing
        elseif node1 == node2.parent
            return node1
        else
            node2 = node2.parent
        end
    end
end

function findMRCA(nodes::Vector)
    nodes = copy(nodes)
    node1 = pop!(nodes)
    while length(nodes) > 0
        node2 = pop!(nodes)
        node1 = findMRCA(node1, node2)
    end
    return node1
end

function findMRCA(population)
    return findMRCA(allcells(population))
end

#endregion

allcells(population) = BinaryNode{SimpleTreeCell}[x for x in population.cells if !isnothing(x)]
