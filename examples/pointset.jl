using PerturbedLattices

pts = PointSet(10,2)

move = GaussianMoveModel([0.5 0.; 0. 0.5])

rand!(pts, move)