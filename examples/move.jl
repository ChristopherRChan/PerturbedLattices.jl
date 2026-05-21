using PerturbedLattices

move = GaussianMoveModel([0.5 0.; 0. 0.5])
rand(move)

move2 = UniformMoveModel([[-1.0, 1.0], [-1.0, 1.0]])
rand(move2)

move = GaussianMoveModel(√0.5, 2)
rand(move)
