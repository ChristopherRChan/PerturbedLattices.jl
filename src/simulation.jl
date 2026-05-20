"""
Monte Carlo simulation functions using Metropolis-Hasting algorithm. 
"""
function iterate!(rng::AbstractRNG, pl::PerturbedLatticeModel)
    sqdist = square_distance(pl.pointset)

    # Proposal of new point 
    i, new_point = proposal(rng, pl)
    # before move
    old_loc_en = local_energy(pl.h, i)
    # after move
    move!(pl, i, new_point)
    new_loc_en = local_energy(pl.h, i)

    # Metropolis-Hastings acceptance ratio
    r = (new_loc_en == Inf && old_loc_en != Inf) ? 0.0 : r = exp(-(new_loc_en - old_loc_en))
    if rand(rng) > r
        revert!(pl, i)
    end
end

function Random.rand!(rng::AbstractRNG, pl::PerturbedLatticeModel; NMC::Int=10000)
    for _ in 1:NMC
        iterate!(rng, pl)
    end
end

Random.rand!(pl::PerturbedLatticeModel; NMC::Int=1000) = Random.rand!(Random.default_rng(), pl; NMC=NMC)