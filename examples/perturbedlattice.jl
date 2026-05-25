using PerturbedLattices

h = StraussHamiltonian(β=1.0, ρ=1.3)
move = GaussianMove(σ²=0.5, d=2)

pl = PerturbedLatticeModel(h, move, (20, 2))
dim(pl)
pl.θ
θ(pl)
rand!(pl.pointset, move)
points(pl)

pl[1]
pl[-10,-10]