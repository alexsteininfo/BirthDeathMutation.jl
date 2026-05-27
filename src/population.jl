abstract type AbstractPopulation end

"""
    Population

Holds a single population of cells as a flat vector of `BinaryNode{SimpleTreeCell}` nodes.
`Nothing` entries are tombstones left after cell death so that tree parent pointers remain
valid; they are removed lazily by `allcells`.
"""
mutable struct Population <: AbstractPopulation
    cells::Vector{Union{BinaryNode{SimpleTreeCell}, Nothing}}
    t::Float64
    subclones::Vector{Subclone}
end

Base.size(population::Population) = (length(population),)
Base.length(population::Population) = length(population.cells)

function Base.iterate(population::Population)
    return iterate(population.cells)
end

function Base.iterate(population::Population, state)
    return iterate(population.cells, state)
end

function Base.getindex(population::Population, i)
    return getindex(population.cells, i)
end

function Base.getindex(population::Population, a::Vector{Int64})
    return map(i -> getindex(population, i), a)
end

function Base.setindex!(population::Population, v, i)
    return setindex!(population.cells, v, i)
end

Base.firstindex(population::Population) = firstindex(population.cells)
Base.lastindex(population::Population) = lastindex(population.cells)

function Population(
    cells,
    birthrate,
    deathrate,
    moranrate,
    asymmetricrate
)
    return Population(
        cells,
        0.0,
        Subclone[Subclone(
            1,
            0,
            0.0,
            length(filter(!isnothing, cells)),
            birthrate,
            deathrate,
            moranrate,
            asymmetricrate
        )]
    )
end

function Base.show(io::IO, population::Population)
    @printf(io, "Population: \n    %d cells", length(allcells(population)))
    @printf(io, "\n    %d subclones", length(filter(x -> x.size > 0, population.subclones)))
    @printf(io, " (t = %.2f)", population.t)
end
