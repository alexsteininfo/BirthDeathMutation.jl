"""
    simulate!(population, input, selection, nextID, rng::AbstractRNG=Random.GLOBAL_RNG;
        timefunc=exptime, t0=nothing, tmax=nothing)

Run a simulation defined by `input` and `selection`, starting from `population`.
"""
function simulate! end

function simulate!(
    population::Population,
    input::BranchingMoranInput,
    selection::AbstractSelection,
    nextID::Integer,
    rng::AbstractRNG=Random.GLOBAL_RNG;
    timefunc=exptime,
    t0=nothing,
    tmax=nothing
)
    population, nextID = branchingprocess!(
        population,
        selection,
        input.Nmax,
        input.μ,
        input.mutationdist,
        isnothing(tmax) ? input.tmax : minimum((tmax, input.tmax)),
        nextID,
        rng;
        timefunc,
        t0
    )
    population, nextID = moranprocess!(
        population,
        selection,
        input.μ,
        input.mutationdist,
        isnothing(tmax) ? input.tmax : minimum((tmax, input.tmax)),
        nextID,
        rng;
        timefunc,
        t0,
        moranincludeself=input.moranincludeself
    )
    final_timedep_mutations!(population, input.μ, input.mutationdist, rng; tend=input.tmax)
    return population, nextID
end

function simulate!(
    population::Population,
    input::BranchingInput,
    selection::AbstractSelection,
    nextID::Integer,
    rng::AbstractRNG=Random.GLOBAL_RNG;
    timefunc=exptime,
    t0=nothing,
    tmax=nothing
)
    population, nextID = branchingprocess!(
        population,
        selection,
        input.Nmax,
        input.μ,
        input.mutationdist,
        isnothing(tmax) ? input.tmax : minimum((tmax, input.tmax)),
        nextID,
        rng;
        timefunc,
        t0
    )
    final_timedep_mutations!(population, input.μ, input.mutationdist, rng; tend=input.tmax)
    return population, nextID
end

function simulate!(
    population::Population,
    input::MoranInput,
    selection::AbstractSelection,
    nextID::Integer,
    rng::AbstractRNG=Random.GLOBAL_RNG;
    timefunc=exptime,
    t0=nothing,
    tmax=nothing
)
    population, nextID = moranprocess!(
        population,
        selection,
        input.μ,
        input.mutationdist,
        isnothing(tmax) ? input.tmax : minimum((tmax, input.tmax)),
        nextID,
        rng;
        timefunc,
        t0,
        moranincludeself=input.moranincludeself
    )
    final_timedep_mutations!(population, input.μ, input.mutationdist, rng; tend=input.tmax)
    return population, nextID
end

"""
    branchingprocess!(population, selection, Nmax, μ, mutationdist, tmax, nextID, rng;
        timefunc=exptime, t0=nothing)

Run a stochastic branching process on `population` until it reaches `Nmax` cells or `tmax`.
"""
function branchingprocess!(
    population::Population,
    selection::AbstractSelection,
    Nmax,
    μ,
    mutationdist,
    tmax,
    nextID,
    rng::AbstractRNG;
    timefunc=exptime,
    t0=nothing
)
    t = !isnothing(t0) ? t0 : age(population)
    N = length(allcells(population))

    nsubclones = getmaxsubclones(selection)
    nsubclonescurrent = length(population.subclones)
    birthrates = getbirthrates(population.subclones)
    deathrates = getdeathrates(population.subclones)
    Rmax = maximum(birthrates) + maximum(deathrates)

    while N < Nmax && N > 0
        Δt = 1 / (Rmax * N) .* timefunc(rng)
        t + Δt <= tmax || break
        t += Δt

        population, birthrates, deathrates, Rmax, N, nextID, nsubclonescurrent, nsubclones =
            branchingupdate!(
                population, selection, birthrates, deathrates, Rmax, N, nextID,
                nsubclonescurrent, nsubclones, t, μ, mutationdist, rng
            )
    end
    return population, nextID
end

function branchingupdate!(
    population::Population,
    selection,
    birthrates,
    deathrates,
    Rmax,
    N,
    nextID,
    nsubclonescurrent,
    nsubclones,
    t,
    μ,
    mutationdist,
    rng
)
    randcellid = rand(rng, 1:N)
    randcell = population.cells[randcellid]
    r = rand(rng, Uniform(0, Rmax))
    cellsubclone = getclonetype(randcell)
    br = birthrates[cellsubclone]
    dr = deathrates[cellsubclone]

    if r < br
        population, _, nextID = celldivision!(
            population, population.subclones, randcellid, t, nextID, μ, mutationdist, rng
        )
        N += 1
        if newsubclone_ready(selection, nsubclonescurrent, nsubclones, t, rng)
            newmutant_selectioncoeff =
                getselectioncoefficient(selection, nsubclonescurrent, rng)
            cellmutation!(
                population, population.subclones, newmutant_selectioncoeff,
                population.cells[randcellid], t
            )
            nsubclonescurrent += 1
            birthrates = getbirthrates(population.subclones)
            deathrates = getdeathrates(population.subclones)
            Rmax = maximum(birthrates) + maximum(deathrates)
        end
    elseif r < br + dr
        population, _ = celldeath!(population, population.subclones, randcellid, t)
        N -= 1
    end
    updatetime!(population, t)
    return population, birthrates, deathrates, Rmax, N, nextID, nsubclonescurrent, nsubclones
end

"""
    moranprocess!(population, selection, μ, mutationdist, tmax, nextID, rng;
        timefunc=exptime, t0=nothing, moranincludeself=true)

Run a Moran process on `population` until `tmax`.
"""
function moranprocess!(
    population::Population,
    selection,
    μ,
    mutationdist,
    tmax,
    nextID,
    rng::AbstractRNG;
    timefunc=exptime,
    t0=nothing,
    moranincludeself=true
)
    t = !isnothing(t0) ? t0 : age(population)
    N = length(allcells(population))
    nsubclones = getmaxsubclones(selection)
    nsubclonescurrent = length(population.subclones)
    moranrates = getmoranrates(population.subclones)
    Rmax = maximum(moranrates)
    while true
        Δt = 1 / (Rmax * N) .* timefunc(rng)
        t = t + Δt
        if t > tmax
            break
        end
        population, moranrates, Rmax, N, nextID, nsubclonescurrent = moranupdate!(
            population, selection, moranrates, Rmax, N, nextID, nsubclonescurrent,
            nsubclones, t, μ, mutationdist, moranincludeself, rng
        )
    end
    return population, nextID
end

function moranupdate!(
    population::Population,
    selection,
    moranrates,
    Rmax,
    N,
    nextID,
    nsubclonescurrent,
    nsubclones,
    t,
    μ,
    mutationdist,
    moranincludeself,
    rng
)
    dividecellid = rand(rng, 1:N)
    dividecell = population.cells[dividecellid]
    r = rand(rng, Uniform(0, Rmax))
    mr = moranrates[getclonetype(dividecell)]
    if r < mr
        deadcell = choose_moran_deadcell(N, dividecellid, moranincludeself, rng)
        population, _, nextID = celldivision!(
            population, population.subclones, dividecellid, t, nextID, μ, mutationdist, rng
        )
        if newsubclone_ready(selection, nsubclonescurrent, nsubclones, t, rng)
            newmutant_selectioncoeff =
                getselectioncoefficient(selection, nsubclonescurrent, rng)
            cellmutation!(
                population, population.subclones, newmutant_selectioncoeff,
                population.cells[dividecellid], t
            )
            nsubclonescurrent += 1
            moranrates = getmoranrates(population.subclones)
            Rmax = maximum(moranrates)
        end
        population, _ = celldeath!(population, population.subclones, deadcell, t)
    end
    updatetime!(population, t)
    return population, moranrates, Rmax, N, nextID, nsubclonescurrent
end

function choose_moran_deadcell(modulesize, dividecellid, moranincludeself, rng)
    if moranincludeself
        deadcellid = rand(rng, 1:modulesize)
        if deadcellid == dividecellid
            return modulesize + 1
        else
            return deadcellid
        end
    else
        return rand(rng, deleteat!(collect(1:modulesize), dividecellid))
    end
end

updatetime!(population::Population, t) = (population.t = t)

function discretetime(rng, λ=1)
    return 1/λ
end

function exptime(rng::AbstractRNG)
    rand(rng, Exponential(1))
end

function exptime(rng::AbstractRNG, λ)
    rand(rng, Exponential(1/λ))
end

age(population::Population) = population.t
age(simulation::Simulation) = age(simulation.output)
