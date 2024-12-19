
using PerturbedLattices
gpl = GibbsianPerturbedLattice(radius=(10.0,10.0), interactions = [pairwise(l -> 4 * ((1 / l)^12 - (1 / l)^6))])
gpl.interactions[1].φ(2.0)
rand(gpl, :inside)
using GeoStats
using CairoMakie 
viz(gpl.points)
length(gpl.points)
gpl