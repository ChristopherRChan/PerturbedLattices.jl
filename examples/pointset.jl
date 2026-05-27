using PerturbedLattices

ps = PointSet(10,2)
dim(ps)
convert(PointSet, ps)

move = GaussianMove([0.5 0.; 0. 0.5])

rand!(ps, move)