"""
    simulate!(population, block::BirthDeathBlock, rng) -> Population
    simulate!(population, block::MoranBlock,       rng) -> Population

Run a simulation block on `population`, extending its tree in-place. Returns the
same `population` so blocks can be chained:

```julia
pop = simulate!(pop, block1, rng)
pop = simulate!(pop, block2, rng)
```
"""
function simulate! end

function simulate!(
    population::Population,
    block::BirthDeathBlock,
    rng::AbstractRNG=Random.GLOBAL_RNG
)
    if block.restart_on_extinction
        initial_cells     = deepcopy(population.cells)
        initial_t         = population.t
        initial_subclones = deepcopy(population.subclones)
    end

    while true
        t                 = age(population)
        N                 = length(allcells(population))
        nextID            = getnextID(population)
        nsubclones        = getmaxsubclones(block.selection)
        nsubclonescurrent = length(population.subclones)

        while !block.stopfunction(population) && N > 0
            Rmax = maximum(
                block.birthrate(sc.s, N) + block.deathrate(sc.s, N)
                for sc in population.subclones
            )
            Δt = exptime(rng) / (Rmax * N)
            t += Δt
            population, N, nextID, nsubclonescurrent = _branchingupdate!(
                population, block, Rmax, N, nextID, nsubclonescurrent, nsubclones, t, rng
            )
        end

        if N == 0 && block.restart_on_extinction
            population.cells     = deepcopy(initial_cells)
            population.t         = initial_t
            population.subclones = deepcopy(initial_subclones)
        else
            break
        end
    end
    return population
end

function simulate!(
    population::Population,
    block::MoranBlock,
    rng::AbstractRNG=Random.GLOBAL_RNG
)
    N = length(allcells(population))
    N == block.N || error(
        "Population size $N does not match MoranBlock.N = $(block.N)"
    )
    t = age(population)
    nextID = getnextID(population)
    nsubclones = getmaxsubclones(block.selection)
    nsubclonescurrent = length(population.subclones)

    while !block.stopfunction(population)
        Rmax = maximum(block.moranrate(sc.s, N) for sc in population.subclones)
        Δt = exptime(rng) / (Rmax * N)
        t += Δt

        population, nextID, nsubclonescurrent = _moranupdate!(
            population, block, Rmax, N, nextID, nsubclonescurrent, nsubclones, t, rng
        )
    end
    return population
end

function _branchingupdate!(
    population, block, Rmax, N, nextID, nsubclonescurrent, nsubclones, t, rng
)
    randcellid = rand(rng, 1:N)
    randcell = population.cells[randcellid]
    cellsubclone = population.subclones[getclonetype(randcell)]
    r = rand(rng, Uniform(0, Rmax))
    br = block.birthrate(cellsubclone.s, N)
    dr = block.deathrate(cellsubclone.s, N)

    if r < br
        population, _, nextID = celldivision!(
            population, population.subclones, randcellid, t, nextID,
            block.μ, block.mutationdist, rng
        )
        N += 1
        if newsubclone_ready(block.selection, nsubclonescurrent, nsubclones, t, rng)
            s = getselectioncoefficient(block.selection, nsubclonescurrent, rng)
            cellmutation!(
                population, population.subclones, s, population.cells[randcellid], t
            )
            nsubclonescurrent += 1
        end
    elseif r < br + dr
        population, _ = celldeath!(population, population.subclones, randcellid, t)
        N -= 1
    end
    updatetime!(population, t)
    return population, N, nextID, nsubclonescurrent
end

function _moranupdate!(
    population, block, Rmax, N, nextID, nsubclonescurrent, nsubclones, t, rng
)
    dividecellid = rand(rng, 1:N)
    dividecell = population.cells[dividecellid]
    r = rand(rng, Uniform(0, Rmax))
    mr = block.moranrate(population.subclones[getclonetype(dividecell)].s, N)

    if r < mr
        deadcellid = _choose_moran_deadcell(N, dividecellid, block.moranincludeself, rng)
        population, _, nextID = celldivision!(
            population, population.subclones, dividecellid, t, nextID,
            block.μ, block.mutationdist, rng
        )
        if newsubclone_ready(block.selection, nsubclonescurrent, nsubclones, t, rng)
            s = getselectioncoefficient(block.selection, nsubclonescurrent, rng)
            cellmutation!(
                population, population.subclones, s, population.cells[dividecellid], t
            )
            nsubclonescurrent += 1
        end
        population, _ = celldeath!(population, population.subclones, deadcellid, t)
    end
    updatetime!(population, t)
    return population, nextID, nsubclonescurrent
end

function _choose_moran_deadcell(N, dividecellid, moranincludeself, rng)
    if moranincludeself
        deadcellid = rand(rng, 1:N)
        return deadcellid == dividecellid ? N + 1 : deadcellid
    else
        return rand(rng, deleteat!(collect(1:N), dividecellid))
    end
end

updatetime!(population::Population, t) = (population.t = t)

function exptime(rng::AbstractRNG)
    rand(rng, Exponential(1))
end

age(population::Population) = population.t
