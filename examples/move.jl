using PerturbedLattices

move = GaussianMove([0.5 0.; 0. 0.5])
rand(move)

move2 = UniformMove([[-1.0, 1.0], [-1.0, 1.0]])
rand(move2)

move = GaussianMove(√0.5, 2)
rand(move)
