using WGLMakie 
using PerturbedLattices
pl = GibbsPerturbedLattice()

function φ(r::Float64) 
if r <= 1
    return 10
else 
    return 0
end
end
H = [PerturbedLattices.pairwise(φ)] 
gpl = GibbsPerturbedLattice(interactions = H, nsim=10^4)
ps = rand(gpl, 10.0)
viz(ps)
