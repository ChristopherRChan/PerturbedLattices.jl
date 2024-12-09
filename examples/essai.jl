using PerturbedLattices
pl = PerturbedLattice()
ps = rand(pl,10.0)
using CairoMakie 
viz(ps)