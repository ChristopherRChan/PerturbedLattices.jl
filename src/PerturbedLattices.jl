module PerturbedLattices

using LinearAlgebra
using Random
using Distributions
using Plots
using Optim
using OffsetArrays

export PerturbedLatticeModel
export PointSet, Grid
export MultiStraussHamiltonian, LennardJonesHamiltonian, StraussHamiltonian, HardCoreHamiltonian
export GaussianMove, UniformMove
export lattice, points, points!
export TakacsFiksel
export rand!
export plot 
export θ, nθ, θ!
export dim
export hγ, d²
export fit!

# Core types
export PerturbedLatticeV1

# Main functions
export create_grid, iterate!, shift!
export local_energy
export points_in_window 
export DLR_W, fit
export plot_points, plot_point_grid_connection

const Point = Vector{Float64}
const Points = Vector{Point}
# const OffsetPoints = OffsetArray{Point}

abstract type AbstractPoints end # has a method points  and a method index for possible negative coordinates
Base.getindex(p::AbstractPoints, i) = Base.getindex(points(p), i)
Base.getindex(p::AbstractPoints, i, j...) = Base.getindex(points(p), index(p)[i, j...])
Base.setindex!(p::AbstractPoints, v, i) = Base.setindex!(points(p), v, i)
Base.setindex!(p::AbstractPoints, v, i, j...) = Base.setindex!(points(p), v, index(p)[i, j...])

abstract type AbstractPointSet <: AbstractPoints end # has a field points
points(ps::AbstractPointSet) = ps.points
points!(ps::AbstractPointSet, points::Points) = ps.points = points
points!(ps::AbstractPointSet, ps2::AbstractPoints) = points!(ps, points(ps2))

abstract type AbstractLattice <: AbstractPointSet end

abstract type AbstractPerturbedLatticeModel <: AbstractPointSet end
abstract type EstimationMethod end


abstract type AbstractHamiltonian <: AbstractPointSet end
abstract type ExponentialFamilyHamiltonian <: AbstractHamiltonian end
abstract type AbstractMove end
points(h::AbstractHamiltonian) = points(h.pointset)

# Include submodules
include("dim.jl")
include("param.jl")
include("grid.jl")
include("pointset.jl")
include("hamiltonian.jl")
include("move.jl")
#include("estimation.jl")
#include("data_creation.jl")
include("perturbed_lattice_model.jl")
include("simulation.jl")
include("visualization.jl")
include("estim_takacs_fiksel.jl")

include("v1/perturbedLatticeModel.jl")
include("v1/grid.jl")
include("v1/energy.jl")
include("v1/simulation.jl")
include("v1/visualization.jl")
include("v1/estimation.jl")
include("v1/data_creation.jl")

end