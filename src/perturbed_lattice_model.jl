"""
    PerturbedLattice

A structure for managing perturbed lattice point processes in 2D or 3D.

# Fields

"""

mutable struct PerturbedLatticeModel <: AbstractPerturbedLatticeModel
    h::AbstractHamiltonian
    move::AbstractMoveModel
    pointset::PointSet
    # internally
    points::Points
    ind::OffsetArray{Int}

    function PerturbedLatticeModel(h::AbstractHamiltonian, move::AbstractMoveModel, pointset::PointSet)
        pl = new(h, move, pointset)
        ## delegate
        pl.points = pl.pointset.points
        pl.ind = pl.pointset.ind
        return pl
    end
end

function PerturbedLatticeModel(h::AbstractHamiltonian, move::AbstractMoveModel, radius::Int = 20, d::Int=2)
    ps = PointSet(radius, d)
    pointset!(h, ps)
    pl = PerturbedLatticeModel(h, move, ps)
    return pl
end

Base.length(pl::PerturbedLatticeModel) = length(pl.pointset.lattice)

function Base.show(io::IO, pl::PerturbedLatticeModel)
    print(io, "PerturbedLatticeModel  on " * string(pl.pointset.lattice))
end

points(pl::PerturbedLatticeModel) = points(pl.pointset)

lattice(pl::PerturbedLatticeModel) = lattice(pl.pointset)


update!(pl::PerturbedLatticeModel) = update!(pl.pointset)
move!(pl::PerturbedLatticeModel, i::Int, point::Point) = move!(pl.pointset, i, point)
revert!(pl::PerturbedLatticeModel, i::Int) = revert!(pl.pointset, i)

function proposal(rng::AbstractRNG, pl::PerturbedLatticeModel; radius::Int=pl.pointset.lattice.radius)
    ii = rand(rng, -radius:radius, pl.pointset.lattice.d)
    i = pl.ind[ii...]
    return (i,lattice(pl)[i] .+ rand(rng, pl.move))
end
   