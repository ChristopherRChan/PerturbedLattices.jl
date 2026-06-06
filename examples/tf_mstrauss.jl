using PerturbedLattices
# using BenchmarkTools


h = MultiStraussHamiltonian(β=[-1.0,2.0], ρ=[0.0,1.0,2.0])
move = GaussianMove(σ²=1.0, d=2)

# Create the lattice
pl = PerturbedLatticeModel(h, move, (20, 2))

# create TF estimation
tf = TakacsFiksel(pl, 15, nQ=10, ρQ=3.0)

## Repeat these 4 lines to have new estimation
θ!(pl, [1.,-1.0,2.0])
@time rand!(pl,NMC=100000)
fit!(tf, [1.0, .2, .2])
tf