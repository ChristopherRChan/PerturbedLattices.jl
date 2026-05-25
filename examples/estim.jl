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
PerturbedLattices.cache!(tf)
points(tf.gridQ)

tf.fcache
tf.ΣfΛcache
tf.ΣfΛcache[lattice(tf.pl).ind[841], : , 2]

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