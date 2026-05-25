"""
    PerturbedLattice

A structure for managing perturbed lattice point processes in 2D or 3D.

# Fields

"""

mutable struct PerturbedLatticeModel <: AbstractPerturbedLatticeModel
    h::AbstractHamiltonian
    move::AbstractMove
    pointset::PointSet
    # internally
    θ::Vector{Float64}      # to keep the current parameters (after estimation) 
    points::Points          # a pointer to the pointset.points
    ind::OffsetArray{Int}   # OffsetArray to convert cartesian with negative index to index of point

    function PerturbedLatticeModel(h::AbstractHamiltonian, move::AbstractMove, pointset::PointSet)
        pl = new(h, move, pointset)
        pl.θ = θ(pl) # at first save parameters of move and hamiltonian  inside internal vector θ
        ## delegate
        pl.points = pl.pointset.points
        pl.ind = pl.pointset.ind
        return pl
    end
end

function PerturbedLatticeModel(h::AbstractHamiltonian, move::AbstractMove, radius::Int = 20, d::Int=2)
    ps = PointSet(radius, d)
    pointset!(h, ps)
    pl = PerturbedLatticeModel(h, move, ps)
    return pl
end

PerturbedLatticeModel(h::AbstractHamiltonian, move::AbstractMove, grid::Tuple{Int64, Int64}) = PerturbedLatticeModel(h, move, grid[1], grid[2])

nbparam(pl::PerturbedLatticeModel) = nbparam(pl.move) + nbparam(pl.h)
params(pl::PerturbedLatticeModel) = [params(pl.move); params(pl.h)]
function params!(pl, θ::Vector{Float64})
    if nbparam(pl.move) > 0
        params!(pl.move, θ[1:nbparam(pl.move)])
    end
    if nbparam(pl.h) > 0
        params!(pl.h, θ[1 + nbparam(pl.move):nbparam(pl)])
    end
end

Base.length(pl::PerturbedLatticeModel) = length(pl.pointset.lattice)

function Base.show(io::IO, pl::PerturbedLatticeModel)
    print(io, "PerturbedLatticeModel  on " * string(pl.pointset.lattice))
end

points(pl::PerturbedLatticeModel) = points(pl.pointset)
points!(pl::PerturbedLatticeModel, points::Points) = pl.points = points
points!(pl::PerturbedLatticeModel, pointset::PointSet) = begin pl.pointset = pointset; pl.points = pointset.points; end

lattice(pl::PerturbedLatticeModel) = lattice(pl.pointset)


update!(pl::PerturbedLatticeModel) = update!(pl.pointset) 

function move!(pl::AbstractPerturbedLatticeModel, i::Int, point::Point)
    # before move
    old_en = hγ(pl.h, i)
    # after move
    move!(pl.pointset, i, point)
    new_en = hγ(pl.h, i)
    return (new_en == Inf && old_en != Inf) ? 0.0 : exp(-(new_en - old_en))
end

revert!(pl::PerturbedLatticeModel, i::Int) = revert!(pl.pointset, i)

function proposal(rng::AbstractRNG, pl::PerturbedLatticeModel; radius::Int=pl.pointset.lattice.radius)
    ii = rand(rng, -radius:radius, pl.pointset.lattice.d)
    i = pl.ind[ii...]
    return (i,lattice(pl)[i] .+ rand(rng, pl.move))
end

# function local_energy(pl::PerturbedLatticeModel)
   