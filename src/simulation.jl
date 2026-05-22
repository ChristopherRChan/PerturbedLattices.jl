"""
Monte Carlo simulation functions using Metropolis-Hasting algorithm. 
"""
function iterate!(rng::AbstractRNG, pl::PerturbedLatticeModel, radius::Int=pl.pointset.lattice.radius)    
    # Proposal of new point 
    i, new_point = proposal(rng, pl; radius=radius)
    # Metropolis-Hastings acceptance ratio
    if rand(rng) > move!(pl, i, new_point)
        revert!(pl, i)
    end
end

function Random.rand!(rng::AbstractRNG, pl::PerturbedLatticeModel, radius::Int=pl.pointset.lattice.radius; NMC::Int=10000 )
    for _ in 1:NMC
        iterate!(rng, pl, radius)
    end
end

Random.rand!(pl::PerturbedLatticeModel, radius::Int=pl.pointset.lattice.radius; NMC::Int=10000) = Random.rand!(Random.default_rng(), pl, radius; NMC=NMC)