mutable struct GibbsianPerturbedLattice <: GeoStatsProcesses.PointProcess
    δ::Float64
    q::Distribution
    interactions::Vector{Interaction}
    nsim::Int64
    radius::Tuple{Vararg{Float64}}
    grid::Vector
    points::PointSet
    
    function GibbsianPerturbedLattice(;δ::Float64=1.0, q::Distribution=Normal(0.1), interactions = [], radius::Tuple{Vararg{Float64}} = (10.0,10.0), nsim::Int64=10^4)
        gpl = new()
        gpl.δ, gpl.q, gpl.nsim = δ, q, nsim
        gpl.interactions = convert.(Interaction, interactions)
        gpl.radius = radius 
        initpoints!(gpl) 
        return gpl
    end
end

function initpoints!(gpl::GibbsianPerturbedLattice)
    d = length(gpl.radius)
    gpl.grid = vertices(CartesianGrid(.-(gpl.radius), gpl.radius; dims = Int.(2 .* gpl.radius ./ gpl.δ)))
    gpl.points = PointSet([ v + Vec(rand(gpl.q, d)...) for v in gpl.grid ])
end


function GeoStatsProcesses.randsingle(rng::AbstractRNG, gpl::GibbsianPerturbedLattice, g) #g=radius of box
    cpt = 0
    d = length(gpl.radius)
    n = length(gpl.points)
    for _ in 1:gpl.nsim
        i = rand(rng, 1:n)
        newpt = gpl.grid[i] + Vec(rand(rng, gpl.q, d)...)
        if rand() <= exp( -moveenergy(gpl, i, newpt) )
            cpt += 1
            gpl.points.geoms[i] = newpt
        end
    end
    println(cpt, " moves over ", gpl.nsim)
    return gpl.points
end

function moveenergy(gpl::GibbsianPerturbedLattice, i::Int64, pt::Point)
    en = 0.0
    for interaction in gpl.interactions
        en += moveenergy(interaction, i, pt, gpl.points)
    end
    return en
end