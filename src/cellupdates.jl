"""
    numbernewmutations(rng, mutationdist, μ; Δt=nothing)

Generate the number of new mutations for a given mutational process `mutationdist`
    with mean `μ` (non-time independent, e.g. :fixed, :poisson) or `μΔt` (time dependent,
    e.g. :poissontimedep).
"""
function numbernewmutations(rng, mutationdist, μ; Δt=nothing)
    if mutationdist == :fixed
        return round(Int64, μ)
    elseif mutationdist == :poisson
        return rand(rng, Poisson(μ))
    elseif mutationdist == :geometric
        return rand(rng, Geometric(1/(1+μ)))
    elseif mutationdist == :poissontimedep
        return rand(rng, Poisson(μ*Δt))
    elseif mutationdist == :fixedtimedep
        return round(Int64, μ*Δt)
    else
        error("$mutationdist is not a valid mutation rule")
    end
end

"""
    celldivision!(population, subclones, parentcellid, t, nextID, μ, mutationdist, rng;
        nchildcells=2)

Cell at `population.cells[parentcellid]` divides. If `nchildcells == 1` it is replaced by a
single child cell. If `nchildcells == 2` a second child cell is appended to the end of
`population.cells`. Cell frequencies are updated in `subclones`.

Time-dependent mutations are assigned to the parent node before division; non-time-dependent
mutations are assigned to the child cells at division.
"""
function celldivision!(
    population::Population,
    subclones,
    parentcellid,
    t,
    nextID,
    μ,
    mutationdist,
    rng;
    nchildcells=2,
    timedepmutationsonly=false
)
    alivecells = population.cells
    parentcellnode = alivecells[parentcellid]
    childcellmuts = zeros(Int64, nchildcells)
    for (μ0, mutationdist0) in zip(μ, mutationdist)
        if mutationdist0 == :fixedtimedep || mutationdist0 == :poissontimedep
            Δt = t - parentcellnode.data.latestupdatetime
            parentcellnode.data.mutations += numbernewmutations(
                rng, mutationdist0, μ0, Δt=Δt
            )
            parentcellnode.data.latestupdatetime = t
        elseif !timedepmutationsonly
            for i in 1:nchildcells
                childcellmuts[i] += numbernewmutations(rng, mutationdist0, μ0)
            end
        end
    end
    childcell1 = SimpleTreeCell(
        id=nextID,
        birthtime=t,
        mutations=childcellmuts[1],
        clonetype=parentcellnode.data.clonetype
    )
    alivecells[parentcellid] = leftchild!(parentcellnode, childcell1)
    if nchildcells == 2
        childcell2 = SimpleTreeCell(
            id=nextID + 1,
            birthtime=t,
            mutations=childcellmuts[2],
            clonetype=parentcellnode.data.clonetype
        )
        push!(alivecells, rightchild!(parentcellnode, childcell2))
        subclones[parentcellnode.data.clonetype].size += 1
    end
    return population, subclones, nextID + nchildcells
end

"""
    cellmutation!(population, subclones, selectioncoefficient, mutatingcell, t)

Cell mutation occurs in `mutatingcell` to produce a new non-neutral subclone with fitness
equal to `f = 1 + selectioncoefficient`. Birth, moran and asymmetric rates are increased
by factor `f` from wild-type; death rate is unchanged.
"""
function cellmutation!(population, subclones, selectioncoefficient, mutatingcell, t)
    subcloneid = length(subclones) + 1
    parentid = getclonetype(mutatingcell)
    wildtype_rates = getwildtyperates(subclones)
    birthrate, deathrate, moranrate, asymmetricrate = get_newsubclone_rates(
        wildtype_rates, selectioncoefficient
    )
    newsubclone = Subclone(
        subcloneid, parentid, t, 1, birthrate, deathrate, moranrate, asymmetricrate
    )
    push!(subclones, newsubclone)
    setclonetype!(mutatingcell, subcloneid)
    if parentid != 0
        subclones[parentid].size -= 1
    end
    return population, subclones
end

"""
    get_newsubclone_rates(wildtype, selectioncoefficient)

Compute new birth, death, moran and asymmetric rates for a new subclone. All `wildtype`
rates are increased by a factor of `(1 + selectioncoefficient)` except death rate.
"""
function get_newsubclone_rates(wildtype, selectioncoefficient)
    return (
        birthrate = wildtype.birthrate * (1 + selectioncoefficient),
        deathrate = wildtype.deathrate,
        moranrate = wildtype.moranrate * (1 + selectioncoefficient),
        asymmetricrate = wildtype.asymmetricrate * (1 + selectioncoefficient)
    )
end

"""
    celldeath!(population, subclones, deadcellid, t)

Cell at `population.cells[deadcellid]` dies and is removed. Cell frequencies are updated
in `subclones`.
"""
function celldeath!(
    population::Population,
    subclones::Vector{Subclone},
    deadcellid,
    t,
    μ=nothing,
    mutationdist=nothing,
    rng=nothing
)
    alivecells = population.cells
    deadcellclonetype = alivecells[deadcellid].data.clonetype
    killcell!(alivecells, deadcellid, t, μ, mutationdist, rng)
    deleteat!(alivecells, deadcellid)
    subclones[deadcellclonetype].size -= 1
    return population, subclones
end

"""
    getnextID(population)
    getnextID(cells)

Get the next usable cell ID.
"""
function getnextID(population::Population)
    return getnextID(allcells(population))
end

function getnextID(cells::Vector{BinaryNode{SimpleTreeCell}})
    nextID = 1
    for cellnode in cells
        if isnothing(cellnode) continue end
        if id(cellnode) + 1 > nextID
            nextID = id(cellnode) + 1
        end
    end
    return nextID
end

getbirthrates(subclones) = Float64[subclone.birthrate for subclone in subclones]
getdeathrates(subclones) = Float64[subclone.deathrate for subclone in subclones]
getmoranrates(subclones) = Float64[subclone.moranrate for subclone in subclones]
getasymmetricrates(subclones) = Float64[subclone.asymmetricrate for subclone in subclones]

function getwildtyperates(population::AbstractPopulation)
    return getwildtyperates(population.subclones)
end

function getwildtyperates(subclones::Vector{Subclone})
    return (
        birthrate = subclones[1].birthrate,
        deathrate = subclones[1].deathrate,
        moranrate = subclones[1].moranrate,
        asymmetricrate = subclones[1].asymmetricrate
    )
end

initialize_counters(population::Population) = getnextID(population)
