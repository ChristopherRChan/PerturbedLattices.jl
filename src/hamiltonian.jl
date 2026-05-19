abstract type AbstractHamiltonian end

dim(h::AbstractHamiltonian) = dim(h.pointset)

mutable struct StraussHamiltonian <: AbstractHamiltonian
    beta::Float64
    radius::Float64 
    # field required for AbstractHamiltonian
    pointset::AbstractPointSet

    StraussHamiltonian(beta::Float64, radius::Float64) = new(beta, radius)
end

function StraussHamiltonian(beta::Float64, radius::Float64, ps::AbstractPointSet)
    sh = StraussHamiltonian(beta, radius)
    pointset!(sh, ps)
    return sh
end

mutable struct HardCoreHamiltonian <: AbstractHamiltonian
    radius::Float64
    # fields required for AbstractHamiltonian
    pointset::AbstractPointSet

    HardCoreHamiltonian(radius::Float64) = new(radius)
    
end

function HardCoreHamiltonian(radius::Float64, ps::AbstractPointSet)
    hch = HardCoreHamiltonian(radius)
    pointset!(hch, ps)
    return hch
end

function pointset!(sh::StraussHamiltonian, ps::AbstractPointSet)
    sh.pointset = convert(PointSet, ps)
    update!(sh.pointset)
end

function pointset!(hch::HardCoreHamiltonian, ps::AbstractPointSet)
    hch.pointset = convert(PointSet, ps)
    update!(hch.pointset)
end

# local_energy(h::StraussHamiltonian, i::Int) = h.beta * sum(h.adjacency[i, :])

# local_energy(h::HardCoreHamiltonian, i::Int) = any(h.adjacency[i, :]) ? Inf : 0.0

local_energy(h::StraussHamiltonian, i::Int) = h.beta * sum(square_distance(h.pointset)[i, :] .<= h.radius^2)

local_energy(h::HardCoreHamiltonian, i::Int) = any(h.adjacency[i, :] .<= h.radius^2) ? Inf : 0.0


# function adjacency_matrix!(h::AbstractHamiltonian)
#     points = h.pointset.points
#     n_points = length(points)
#     adjacency = falses(n_points, n_points)

#     for i in 1:n_points
#         for j in (i+1):n_points
#             dist_sq = sum((points[i] .- points[j]).^2)
#             if dist_sq <= h.radius^2
#                 adjacency[i, j] = adjacency[j, i] = true
#             end
#         end
#     end

#     h.adjacency = adjacency
# end


# function adjacency_matrix!(h::AbstractHamiltonian, i::Int, new_point::Point)
#     points = h.pointset.points
#     n_points = length(points)

#     # Réinitialiser la ligne et la colonne i
#     h.adjacency[i, :] .= h.adjacency[:, i] .= false

#     # Recalculer les adjacences pour le point i à sa nouvelle position
#     for j in 1:n_points
#         if j != i
#             dist_sq = sum((new_point .- points[j]).^2)
#             if dist_sq <= h.radius^2
#                 h.adjacency[i, j] = h.adjacency[j, i] = true
#             end
#         end
#     end

# end