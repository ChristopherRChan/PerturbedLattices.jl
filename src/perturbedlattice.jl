
struct PerturbedLattice <: GeoStatsProcesses.PointProcess
    δ::Float64
    q::Distribution
end

PerturbedLattice(;δ::Float64=1.0, q::Distribution=Normal(0.1)) = PerturbedLattice(δ, q)

function GeoStatsProcesses.randsingle(rng::AbstractRNG, p::PerturbedLattice, g::Float64) #g=radius of box
    println(2g)
    return PointSet([ v + Vec(rand(p.q, 2)...) for v in vertices(CartesianGrid(Int.((2g,2g)),(0,0), (p.δ, p.δ))) ])
end