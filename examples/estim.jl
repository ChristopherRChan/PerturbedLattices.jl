using PerturbedLattices
# using BenchmarkTools


h = StraussHamiltonian(β=1.0, ρ=2.0)
move = GaussianMove(σ²=0.5, d=2)

# Create the lattice
pl = PerturbedLatticeModel(h, move, (20, 2))

## Warmup phase
println("Starting warmup ...")
θ!(pl, [2,1.0])
@time rand!(pl,NMC=100000)

tf = TakacsFiksel(pl, 15)
PerturbedLattices.cache!(tf)
fit!(tf, [2.0, 2])
tf
@time PerturbedLattices.contrast(tf)
@time PerturbedLattices.contrast(tf, Vector)


### test for devel

dim(tf)
tf.f
points(tf.gridQ)
points(tf)

θ!(tf, tf.pl.θ)
θ(tf)
PerturbedLattices.contrast(tf)
PerturbedLattices.contrast(tf, tf.pl.θ)
fit!(tf, [1.0, 2.])
PerturbedLattices.contrast(tf)
tf.fcache
tf.ΣfΛcache
tf.ΣfΛcache[1, : , :]
lattice(tf.pl).ind[841]
θ(tf)

tf.fcache - [ sum((exp.(-tf.ΣfΛcache[i, : , :]  * θ(tf)) .* tf.ΣfΛcache[i, : , l])) / sum(exp.(-tf.ΣfΛcache[i, : , :]  * θ(tf))) for i=eachindex(tf.subind), l=eachindex(tf.f)]

[ (transpose(fill(1.0, length(tf.subind))) * exp.(-tf.ΣfΛcache[i, : , :]  * θ(tf)) .* tf.ΣfΛcache[i, : , l]) ./ sum(exp.(-tf.ΣfΛcache[i, : , :]  * θ(tf)))  for i=eachindex(tf.subind), l=eachindex(tf.f) ]



transpose(fill(1.0, length(tf.subind))) * 
    (
        tf.fcache - [ sum((exp.(-tf.ΣfΛcache[i, : , :]  * θ(tf)) .* tf.ΣfΛcache[i, : , l])) / sum(exp.(-tf.ΣfΛcache[i, : , :]  * θ(tf))) for i=eachindex(tf.subind), l=eachindex(tf.f)]
        #tf.fcache -  [ sum(exp.(-tf.ΣfΛcache[i, : , :]  * θ(tf)) .* tf.ΣfΛcache[i, : , l]) ./ sum(exp.(-tf.ΣfΛcache[i, : , :]  * θ(tf)))  for i=eachindex(tf.subind), l=eachindex(tf.f) ]
    ) / length(tf.subind)

sum(
    (
        transpose(fill(1.0, length(tf.subind))) * 
        (
            tf.fcache - 
            vcat( [collect(transpose(fill(1.0, length(tf.gridQ))) * (exp.(-tf.ΣfΛcache[i, : , :]  * θ(tf)) .* tf.ΣfΛcache[i, : , :]) / sum(exp.(-tf.ΣfΛcache[i, : , :]  * θ(tf))))  for i=eachindex(tf.subind) ]...)
        ) / length(tf.subind)
    ) .^2
)

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