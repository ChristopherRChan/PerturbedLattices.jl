module PerturbedLattices

using LinearAlgebra
using Random
using Distributions
using StaticArrays
using Plots
using Optim
using JLD2
using Infinity

export PerturbedLatticeModel
export PointSet, Grid
export StraussHamiltonian, HardCoreHamiltonian
export GaussianMoveModel, UniformMoveModel
export lattice, points
export rand!
export plot 

# Core types
export PerturbedLatticeV1

# Main functions
export create_grid, iterate!, shift!
export local_energy
export points_in_window 
export DLR_W, fit
export plot_points, plot_point_grid_connection

const Point= Vector{Float64}
abstract type AbstractPointSet end # has a field points of type Vector{pointset}
abstract type AbstractLattice <: AbstractPointSet end

points(ps::AbstractPointSet) = ps.points

abstract type AbstractPerturbedLatticeModel end

# Include submodules
include("grid.jl")
include("pointset.jl")
include("hamiltonian.jl")
include("move_model.jl")
#include("estimation.jl")
#include("data_creation.jl")
include("perturbed_lattice_model.jl")
include("simulation.jl")
include("visualization.jl")


include("v1/perturbedLatticeModel.jl")
include("v1/grid.jl")
include("v1/energy.jl")
include("v1/simulation.jl")
include("v1/visualization.jl")
include("v1/estimation.jl")
include("v1/data_creation.jl")

end