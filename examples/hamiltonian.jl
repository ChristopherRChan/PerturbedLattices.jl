using PerturbedLattices

h = StraussHamiltonian(1.0, 1.5, Grid(20, 2))

hc = HardCoreHamiltonian(0.5, Grid(20, 2))

local_energy(h, 1)

# constructor without lattice at creation
h2 = StraussHamiltonian(1.0, 1.5)
# lattice initialized after creation
PerturbedLattices.pointset!(h2, Grid(20,2))