mutable struct GibbsPerturbedLattice <: GeoStatsProcesses.PointProcess
    δ::Float64
    q::Distribution
    interactions::Vector{Interaction}
    nsim::Int64
    radius::Tuple{Vararg{Float64}}
    grid::Vector
    points::PointSet
    
    function GibbsPerturbedLattice(;δ::Float64=1.0, q::Distribution=Normal(0.1), interactions = [], radius::Tuple{Vararg{Float64}} = (10.0,10.0), nsim::Int64=10^4)
        gpl = new()
        gpl.δ, gpl.q, gpl.nsim = δ, q, nsim
        gpl.interactions = convert.(Interaction, interactions)
        gpl.radius = radius 
        initpoints!(gpl) 
        return gpl
    end
end

function initpoints!(gpl::GibbsPerturbedLattice)
    d = length(gpl.radius)
    gpl.grid = vertices(CartesianGrid(.-(gpl.radius), gpl.radius; dims = Int.(2 .* gpl.radius ./ gpl.δ)))
    gpl.points = PointSet([ v + Vec(rand(gpl.q, d)...) for v in gpl.grid ])
end