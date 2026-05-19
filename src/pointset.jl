mutable struct PointSet <: AbstractPointSet
    lattice::AbstractLattice
    points::Vector{Point}
end

function PointSet(radius::Int, d::Int) # radius et d définissent la fenêtre d'observation
    g = Grid(radius, d)
    points = deepcopy(g.points) # initalization instead of Vector{Point}(undef, (2*radius+1)^d)
    return PointSet(g, points)
end

points(ps::PointSet) = ps.points

lattice(ps::PointSet) = ps.lattice

dim(ps::PointSet) = ps.lattice.d