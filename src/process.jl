function final_timedep_mutations!(
    population::Population,
    μ,
    mutationdist,
    rng;
    tend=age(population)
)
    for (μ0, mutationdist0) in zip(μ, mutationdist)
        if mutationdist0 ∈ (:poissontimedep, :fixedtimedep)
            for cell in allcells(population)
                Δt = tend - cell.data.latestupdatetime
                cell.data.mutations += numbernewmutations(rng, mutationdist0, μ0, Δt=Δt)
                cell.data.latestupdatetime = tend
            end
        end
    end
end
