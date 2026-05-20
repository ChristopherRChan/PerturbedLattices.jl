using PerturbedLattices
# using BenchmarkTools


h = StraussHamiltonian(1.0, 1.3)
move = GaussianMoveModel([0.5 0.; 0. 0.5])

# Create the lattice
pl = PerturbedLatticeModel(h, move, 20, 2)

## Warmup phase
println("Starting warmup ...")
@time rand!(pl)

println("Warmup completed!\n")
points(pl)
pl[-10,10]
pl[0,0]
# Plot the point grid connection
p = plot(pl)
display(p)

@time rand!(pl, 10)
plot(pl)


println("All simulations completed!")