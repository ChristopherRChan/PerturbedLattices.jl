
struct PerturbedLattice <: GeoStatsProcesses.PointProcess
    δ::Float64
    q::Distribution
end

PerturbedLattice(;δ::Float64=1.0, q::Distribution=Normal(0.1)) = PerturbedLattice(δ, q)

function GeoStatsProcesses.randsingle(rng::AbstractRNG, p::PerturbedLattice, g) #g=radius of box
    d = length(g)
    return PointSet([ v + Vec(rand(rng, p.q, d)...) for v in vertices(CartesianGrid(.-(g), g; dims = Int.(2 .* g ./ p.δ))) ])
end