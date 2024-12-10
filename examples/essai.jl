using GeoStats
using WGLMakie 
using PerturbedLattices
pl = PerturbedLattice()
ps = rand(pl,10.0)
viz(ps)