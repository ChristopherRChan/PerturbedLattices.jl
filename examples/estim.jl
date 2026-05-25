using PerturbedLattices
# using BenchmarkTools


h = StraussHamiltonian(1.0, 1.3)
move = GaussianMoveModel([0.5 0.; 0. 0.5])

# Create the lattice
pl = PerturbedLatticeModel(h, move, 20, 2)

## Warmup phase
println("Starting warmup ...")
@time rand!(pl)

tf = TakacsFiksel(pl, 15)
PerturbedLattices.init!(tf)
tf.f
PerturbedLattices.cache!(tf)

