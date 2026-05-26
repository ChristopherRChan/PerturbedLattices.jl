using PerturbedLattices
# using BenchmarkTools


h = StraussHamiltonian(β=1.0, ρ=1.3)
move = GaussianMove(σ²=0.5, d=2)

# Create the lattice
pl = PerturbedLatticeModel(h, move, (20, 2))

## Warmup phase
println("Starting warmup ...")
@time rand!(pl)

tf = TakacsFiksel(pl, 15)
dim(tf)
tf.f
points(tf.gridQ)
points(tf)

tf.fcache
tf.ΣfΛcache
tf.ΣfΛcache[1, : , :]
lattice(tf.pl).ind[841]
θ(tf)

# left
tf.fcache
# right
transpose(fill(1.0, length(tf.gridQ))) * exp.(-tf.ΣfΛcache[1, : , :]  * θ(tf))
sum(exp.(-tf.ΣfΛcache[1, : , :]  * θ(tf)))
exp.(-tf.ΣfΛcache[1, : , :]  * θ(tf)) .* tf.ΣfΛcache[1, : , :]
collect(transpose(fill(1.0, length(tf.gridQ))) * (exp.(-tf.ΣfΛcache[1, : , :]  * θ(tf)) .* tf.ΣfΛcache[1, : , :]) / sum(exp.(-tf.ΣfΛcache[1, : , :]  * θ(tf))))


vcat([collect(transpose(fill(1.0, length(tf.gridQ))) * (exp.(-tf.ΣfΛcache[i, : , :]  * θ(tf)) .* tf.ΣfΛcache[i, : , :]) / sum(exp.(-tf.ΣfΛcache[i, : , :]  * θ(tf)))) for i=eachindex(tf.subind)]...)

sum(transpose(fill(1.0, length(tf.subind))) * (tf.fcache - vcat([collect(transpose(fill(1.0, length(tf.gridQ))) * (exp.(-tf.ΣfΛcache[i, : , :]  * θ(tf)) .* tf.ΣfΛcache[i, : , :]) / sum(exp.(-tf.ΣfΛcache[i, : , :]  * θ(tf)))) for i=eachindex(tf.subind)]...)) .^2)

#################################################
### Devel of TakacsFiksel
nθ(tf.pl.move)
nθ(tf.pl.h)


d = tf.pl.pointset.lattice.d
tf.radius
subind = lattice(tf.pl).ind[repeat([-tf.radius:tf.radius],d)...]
lattice(tf.pl).ind[0,0]
lattice(tf.pl)[841]
lattice(tf.pl).ind[-tf.radius,tf.radius]
lattice(tf.pl)[1441]
lpts = points(lattice(tf.pl))
lpts[841]

nQ::Int=5
ρQ::Float64=1.5
gridQ = Grid(nQ, d)
δ = ρQ / nQ
PerturbedLattices.scale!(gridQ, δ)
points(gridQ)
gridQ[0,0]
gridQ[-nQ,-nQ]

eachindex(points(tf.gridQ))
eachindex(tf.f)
k=1
si = 841
gridQ[k]
tf.gridQ[k]
pt = lpts[si] .+ gridQ[k]
tf.f[1](si, pt, tf.pl)
tf.f[2](si, pt, tf.pl)