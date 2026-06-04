"""
    initialize_population(N::Int=1; n_mutations::Int=0, time::Float64=0.0) -> Population

Create a population of `N` identical cells, each carrying `n_mutations` mutations,
all belonging to a single wildtype subclone (s=0).
"""
function initialize_population(N::Int=1; n_mutations::Int=0, time::Float64=0.0)
    cells = Vector{Union{BinaryNode{SimpleTreeCell}, Nothing}}(
        [newcell(id, n_mutations) for id in 1:N]
    )
    return Population(cells; time=time)
end

"""
    initialize_population(clones::Vector; time::Float64=0.0) -> Population

Create a population from a vector of clone specifications. Each element must be a
`NamedTuple` with fields `n` (cell count), `s` (selection coefficient), and
`mutations` (mutations per cell). Example:

```julia
pop = initialize_population([
    (n=900, s=0.0, mutations=50),
    (n=100, s=0.1, mutations=55),
])
```
"""
function initialize_population(clones::Vector; time::Float64=0.0)
    cells = Vector{Union{BinaryNode{SimpleTreeCell}, Nothing}}()
    subclones = Subclone[]
    id = 1
    for (i, clone) in enumerate(clones)
        push!(subclones, Subclone(
            subcloneid=i,
            parentid=0,
            mutationtime=0.0,
            size=clone.n,
            s=Float64(clone.s),
            asymmetricrate=0.0
        ))
        for _ in 1:clone.n
            push!(cells, newcell(id, clone.mutations, i))
            id += 1
        end
    end
    return Population(cells, time, subclones)
end

function newcell(id, mutations, clonetype=1)
    return BinaryNode{SimpleTreeCell}(SimpleTreeCell(; id, mutations, clonetype))
end
