"""
    numbernewmutations(rng, mutationdist, μ)

Generate the number of new mutations for a given `mutationdist` (:fixed, :poisson,
or :geometric) with mean `μ`.
"""
function numbernewmutations(rng, mutationdist::Symbol, μ::Float64)
    if mutationdist == :fixed
        return round(Int64, μ)
    elseif mutationdist == :poisson
        return rand(rng, Poisson(μ))
    elseif mutationdist == :geometric
        return rand(rng, Geometric(1 / (1 + μ)))
    else
        error("$mutationdist is not a valid mutation distribution. Use :fixed, :poisson, or :geometric.")
    end
end

"""
    celldivision!(population, subclones, parentcellid, t, nextID, μ, mutationdist, rng;
        nchildcells=2)

Cell at `population.cells[parentcellid]` divides. If `nchildcells == 1` it is replaced
by a single child; if `nchildcells == 2` a second child is appended. Subclone sizes
are updated accordingly.
"""
function celldivision!(
    population::Population,
    subclones,
    parentcellid,
    t,
    nextID,
    μ::Float64,
    mutationdist::Symbol,
    rng;
    nchildcells=2
)
    alivecells = population.cells
    parentcellnode = alivecells[parentcellid]
    childcell1 = SimpleTreeCell(
        id=nextID,
        birthtime=t,
        mutations=numbernewmutations(rng, mutationdist, μ),
        clonetype=parentcellnode.data.clonetype
    )
    alivecells[parentcellid] = leftchild!(parentcellnode, childcell1)
    if nchildcells == 2
        childcell2 = SimpleTreeCell(
            id=nextID + 1,
            birthtime=t,
            mutations=numbernewmutations(rng, mutationdist, μ),
            clonetype=parentcellnode.data.clonetype
        )
        push!(alivecells, rightchild!(parentcellnode, childcell2))
        subclones[parentcellnode.data.clonetype].size += 1
    end
    return population, subclones, nextID + nchildcells
end

"""
    cellmutation!(population, subclones, s, mutatingcell, t)

Cell mutation occurs in `mutatingcell`, creating a new subclone with selection
coefficient `s`. The mutating cell's clonetype is updated to the new subclone.
"""
function cellmutation!(population, subclones, s::Float64, mutatingcell, t)
    subcloneid = length(subclones) + 1
    parentid = getclonetype(mutatingcell)
    parentasymmetricrate = subclones[parentid].asymmetricrate
    newsubclone = Subclone(
        subcloneid=subcloneid,
        parentid=parentid,
        mutationtime=t,
        size=1,
        s=s,
        asymmetricrate=parentasymmetricrate
    )
    push!(subclones, newsubclone)
    setclonetype!(mutatingcell, subcloneid)
    subclones[parentid].size -= 1
    return population, subclones
end

"""
    celldeath!(population, subclones, deadcellid, t)

Cell at `population.cells[deadcellid]` dies and is removed from the population.
"""
function celldeath!(
    population::Population,
    subclones::Vector{Subclone},
    deadcellid,
    t
)
    alivecells = population.cells
    deadcellclonetype = alivecells[deadcellid].data.clonetype
    killcell!(alivecells, deadcellid, t)
    deleteat!(alivecells, deadcellid)
    subclones[deadcellclonetype].size -= 1
    return population, subclones
end

"""
    getnextID(population)

Get the next usable cell ID.
"""
function getnextID(population::Population)
    return getnextID(allcells(population))
end

function getnextID(cells::Vector{BinaryNode{SimpleTreeCell}})
    nextID = 1
    for cellnode in cells
        if !isnothing(cellnode) && id(cellnode) + 1 > nextID
            nextID = id(cellnode) + 1
        end
    end
    return nextID
end

initialize_counters(population::Population) = getnextID(population)
