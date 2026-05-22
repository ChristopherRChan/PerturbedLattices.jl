mutable struct PointSet <: AbstractPointSet
    lattice::AbstractLattice
    points::Points
    ## used internally
    d²::Matrix{Float64}
    old_d²::Vector{Float64}
    old_point::Point
    ind::OffsetArray{Int}

    function PointSet(lattice::AbstractLattice, points::Points)
        ps = new(lattice, points)
        ps.d² = fill(NaN, length(points), length(points))
        ps.ind = ps.lattice.ind
        update!(ps)
        return ps
    end
end

function PointSet(radius::Int, d::Int) # radius et d définissent la fenêtre d'observation
    g = Grid(radius, d)
    points = deepcopy(g.points) # initalization instead of Points(undef, (2*radius+1)^d)
    
    return PointSet(g, points)
end

Base.convert(::Type{PointSet}, g::Grid) = PointSet(g.radius, g.d)

d²(ps::PointSet) = ps.d²
d²(ps::PointSet, i::Int) = ps.d²[i,:]
Σd²(ps::PointSet, i::Int, ρ²::Float64) = sum(d²(ps,i) .<= ρ²)

# iₒ is the index of grid point not considered
d²(ps::PointSet, i::Int, iₒ::Int) = (ps.d²[i,k] for k=eachindex(ps.d²[i,:]) if k ≠ iₒ)
Σd²(ps::PointSet, i::Int, iₒ::Int, ρ²::Float64) = sum(d²(ps,i, iₒ) .<= ρ²)

# Only square distance to update
update!(ps::PointSet) = d²!(ps)

function move!(ps::PointSet, i::Int, point::Point)
    # save the distance of old ith point
    ps.old_d² = copy(ps.d²[i, :])
    ps.old_point = copy(ps.points[i])
    # new point
    ps.points[i] = point
    # update square distance matrix for ith point
    d²!(ps, i)
end

function revert!(ps::PointSet, i::Int)
   ps.points[i] = ps.old_point
   # revert square distance matrix for ith old point
   ps.d²[i, :] = ps.d²[:, i] = ps.old_d²
end

function d²!(ps::PointSet)
    pts = ps.points
    n_pts = length(pts)
    for i in 1:n_pts
        for j in (i+1):n_pts
            @inbounds ps.d²[i, j] = ps.d²[j, i] = sum((pts[i] .- pts[j]).^2)
        end
    end
end

function d²!(ps::PointSet, i::Int)
    pts = ps.points
    n_pts = length(pts)
     for j in 1:n_pts
        if j != i
            @inbounds ps.d²[i, j] = ps.d²[j, i] = sum((pts[i] .- pts[j]).^2)
        end
    end
end

points(ps::PointSet) = ps.points

lattice(ps::PointSet) = ps.lattice

dim(ps::PointSet) = ps.lattice.d