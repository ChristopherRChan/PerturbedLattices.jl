using PerturbedLattices

function φ(r::Float64) 
    if r <= 1
        return 10
    else 
        return 0
    end
end
interact = PerturbedLattices.pairwise(φ)
H = [interact] 
gpl = GibbsPerturbedLattice(interactions = H, nsim=10^4)
ps = rand(gpl, 10.0)
perturbation = ps .- gpl.grid

function ψ(r::Float64)::Float64
    if r <= 1
        return 1
    else 
        return 0
    end
end
interact2 = PerturbedLattices.pairwise(ψ)
H2 = [interact2]
gpl2 = GibbsPerturbedLattice(interactions = H2, nsim=10^4)
fit(gpl2, 1.0, 1.0, ps, perturbation, estimationmethod=true)
