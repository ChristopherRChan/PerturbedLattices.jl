using PerturbedLattices

h = MultiStraussHamiltonian([1.0,-2.0], [0.0,1.0,2.0], Grid(20,2))
hγ(h, 1)

h = LennardJonesHamiltonian(1.0, 1.0, Grid(20,2))
hγ(h, 1, [1.1, 2.1])
PerturbedLattices.S(h, 1, [1.1, 2.1])

h = StraussHamiltonian(2.0, 1.5, Grid(20, 2))

hc = HardCoreHamiltonian(0.5, Grid(20, 2))
h.pointset
d²(h.pointset)
hγ(h, 1)

# constructor without lattice at creation
h2 = StraussHamiltonian(1.0, 1.5)
# lattice initialized after creation
PerturbedLattices.pointset!(h2, Grid(20,2))
h2

dim(h2)