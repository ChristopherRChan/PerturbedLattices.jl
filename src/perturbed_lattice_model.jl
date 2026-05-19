"""
    PerturbedLattice

A structure for managing perturbed lattice point processes in 2D or 3D.

# Fields

"""

mutable struct PerturbedLatticeModel <: AbstractPerturbedLatticeModel
    h::AbstractHamiltonian
    move::AbstractMoveModel
    pointset::PointSet
end

function PerturbedLatticeModel(h::AbstractHamiltonian, move::AbstractMoveModel, radius::Int, d::Int=2)
    ps = PointSet(radius, d)
    lattice!(h, lattice(ps))
    pl = PerturbedLatticeModel(h, move, ps)
    return pl
end

Base.length(pl::PerturbedLatticeModel) = length(pl.pointset.lattice)

function Base.show(io::IO, pl::PerturbedLatticeModel)
    print(io, "PerturbedLatticeModel  on " * string(pl.pointset.lattice))
end

points(pl::PerturbedLatticeModel) = points(pl.pointset)
lattice(pl::PerturbedLatticeModel) = lattice(pl.pointset)
