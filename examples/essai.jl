using GeoStats
using WGLMakie 
using PerturbedLattices
pl = GibbsianPerturbedLattice()
ps = rand(pl,(10.0,10.0,10.0))
viz(ps)