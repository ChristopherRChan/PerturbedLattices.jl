using PerturbedLattices
# using BenchmarkTools


h = LennardJonesHamiltonian(.9, 1.1)
move = GaussianMove(σ²=1.0, d=2)

# Create the lattice
pl = PerturbedLatticeModel(h, move, (20, 2))

# create TF estimation
tf = TakacsFiksel(pl, 15, nQ=5, ρQ=3.0)


## Repeat these 4 lines to have new estimation
pl.θ
θ!(pl, pl.θ)
rand!(pl.pointset, pl.move)
@time rand!(pl,NMC=100000)
plot(pl, radius=3)
fit!(tf, [1.0, 1.0, 1.0])
tf
θ(pl)
θ(move)
θ(h)