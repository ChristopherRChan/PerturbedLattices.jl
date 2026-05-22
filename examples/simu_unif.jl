using PerturbedLattices


h = StraussHamiltonian(-1.0, 1.3)
move = UniformMoveModel([[-0.5, 0.5], [-0.5, 0.5]])

# Create the lattice
pl = PerturbedLatticeModel(h, move, 20, 2)

# Warmup phase
println("Starting warmup ...")
@time rand!(pl, NMC=100000)

println("Warmup completed!\n")
points(pl)
# Plot the point grid connection
p = plot(pl, radius=5)

PerturbedLattices.d²(pl.pointset)

println("All simulations completed!")