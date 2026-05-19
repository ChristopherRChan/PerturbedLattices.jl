using PerturbedLattices

h = StraussHamiltonian(1.0, 1.3)
move = GaussianMoveModel([0.5 0.; 0. 0.5])

pl = PerturbedLatticeModel(h, move, 20)

rand!(pl.pointset, move)
points(pl)