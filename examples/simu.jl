using PerturbedLattices
# using BenchmarkTools


h = StraussHamiltonian(β=1.0, ρ=1.3)
move = GaussianMove(σ²=0.5, d=2)

# Create the lattice
pl = PerturbedLatticeModel(h, move, (20, 2))
θ(pl)
## Warmup phase
println("Starting warmup ...")
@time rand!(pl)
pl.θ
println("Warmup completed!\n")
points(pl)
pl[-10,10]
pl[0,0]
# Plot the point grid connection
plot(pl)

@time rand!(pl, 10)
plot(pl)
plot(pl, radius=10)

println("All simulations completed!")