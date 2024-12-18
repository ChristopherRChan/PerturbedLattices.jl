struct GibbsianPerturbedLattice <: GeoStatsProcesses.PointProcess
    δ::Float64
    q::Distribution
    H::Function
    N::Int64
end

GibbsianPerturbedLattice(;δ::Float64=1.0, q::Distribution=Normal(0.1), H::Function=(pl::Vector -> 0), N::Int64=10^3) = GibbsianPerturbedLattice(δ, q, H, N)

function GeoStatsProcesses.randsingle(rng::AbstractRNG, pl::GibbsianPerturbedLattice, g) #g=radius of box
    d = length(g)
    grid = vertices(CartesianGrid(.-(g), g; dims = Int.(2 .* g ./ pl.δ)))
    ps = [ v + Vec(rand(rng, pl.q, d)...) for v in grid ]
    println(first(ps,3))
    n = length(ps)
    E = 0
    for _ in 1:pl.N
        j=rand(rng, 1:n)
        oldpt =  ps[j]
        newpt = grid[j] + Vec(rand(rng, pl.q, d)...)
        nE = 0
        if rand() <= exp(-nE+E)
            vecps = tempvecPs
            E = nE
        end
    end
    ps = PointSet(vecps)
    return(ps)
end