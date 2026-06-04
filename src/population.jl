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

function Population(cells; time::Float64=0.0, s::Float64=0.0, asymmetricrate::Float64=0.0)
    return Population(
        cells,
        time,
        Subclone[Subclone(
            subcloneid=1,
            parentid=0,
            mutationtime=0.0,
            size=length(filter(!isnothing, cells)),
            s=s,
            asymmetricrate=asymmetricrate
        )]
    )
end

popsize(population::Population) = length(allcells(population))

function Base.show(io::IO, pop::Population)
    @printf(io, "Population: \n    %d cells", length(allcells(pop)))
    @printf(io, "\n    %d subclones", length(filter(x -> x.size > 0, pop.subclones)))
    @printf(io, " (t = %.2f)", pop.t)
end
