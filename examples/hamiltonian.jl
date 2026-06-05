using PerturbedLattices

hMS = MultiStraussHamiltonian([1.0,-2.0], [0.0,1.0,2.0])
hγ(hMS, 1, PointSet(20,2))

hLJ = LennardJonesHamiltonian(1.0, 1.0)
hγ(hLJ, 1, [1.1, 2.1], PointSet(20,2))
PerturbedLattices.S(hLJ, 1, [1.1, 2.1], PointSet(20,2))

hS = StraussHamiltonian(2.0, 1.5)
hγ(hS, 2, PointSet(20,2))
hγ(hS, 2, [1.1,2.1], PointSet(20,2))
PerturbedLattices.S(hS, 2, [1.1,2.1], PointSet(20,2))
hM = MultiStraussHamiltonian([2.0], [0,1.5])
hγ(hM, 2, PointSet(20,2))
hγ(hM, 2, [1.1,2.1], PointSet(20,2))
PerturbedLattices.S(hM, 2, [1.1,2.1], PointSet(20,2))

hc = HardCoreHamiltonian(0.5)
