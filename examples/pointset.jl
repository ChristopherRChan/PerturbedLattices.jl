using PerturbedLattices

ps = PointSet(10,2)
convert(PointSet, ps)

move = GaussianMoveModel([0.5 0.; 0. 0.5])

rand!(ps, move)