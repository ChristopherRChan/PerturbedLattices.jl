using PerturbedLattices

h = StraussHamiltonian(2.0, 1.5, Grid(20, 2))

hc = HardCoreHamiltonian(0.5, Grid(20, 2))
h.pointset
PerturbedLattices.square_distance(h.pointset)
local_energy(h, 1)

# constructor without lattice at creation
h2 = StraussHamiltonian(1.0, 1.5)
# lattice initialized after creation
PerturbedLattices.pointset!(h2, Grid(20,2))
h2