mutable struct PointSet <: AbstractPointSet
    lattice::AbstractLattice
    points::Vector{Point}
    ## used internally
    sqdist::Matrix{Float64}
    old_sqdist::Vector{Float64}
    old_point::Point

    function PointSet(lattice::AbstractLattice, points::Vector{Point})
        ps = new(lattice, points)
        ps.sqdist = fill(NaN, length(points), length(points))
        update!(ps)
        return ps
    end
end

function PointSet(radius::Int, d::Int) # radius et d définissent la fenêtre d'observation
    g = Grid(radius, d)
    points = deepcopy(g.points) # initalization instead of Vector{Point}(undef, (2*radius+1)^d)
    return PointSet(g, points)
end

Base.convert(::Type{PointSet}, g::Grid) = PointSet(g.radius, g.d)

square_distance(ps::PointSet) = ps.sqdist

# Only square distance to update
update!(ps::PointSet) = square_distance!(ps)

function move!(ps::PointSet, i::Int, point::Point)
    # save the distance of old ith point
    ps.old_sqdist = copy(ps.sqdist[i, :])
    ps.old_point = copy(ps.points[i])
    # new point
    ps.points[i] = point
    # update square distance matrix for ith point
    square_distance!(ps, i)
end

function revert!(ps::PointSet, i::Int)
   ps.points[i] = ps.old_point
   # revert square distance matrix for ith old point
   ps.sqdist[i, :] = ps.sqdist[:, i] = ps.old_sqdist
end

function square_distance!(ps::PointSet)
    pts = ps.points
    n_pts = length(pts)
    for i in 1:n_pts
        for j in (i+1):n_pts
            @inbounds ps.sqdist[i, j] = ps.sqdist[j, i] = sum((pts[i] .- pts[j]).^2)
        end
    end
end

function square_distance!(ps::PointSet, i::Int)
    pts = ps.points
    n_pts = length(pts)
     for j in 1:n_pts
        if j != i
            @inbounds ps.sqdist[i, j] = ps.sqdist[j, i] = sum((pts[i] .- pts[j]).^2)
        end
    end
end

points(ps::PointSet) = ps.points

lattice(ps::PointSet) = ps.lattice

dim(ps::PointSet) = ps.lattice.d