abstract type AbstractHamiltonian end

mutable struct StraussHamiltonian <: AbstractHamiltonian
    beta::Float64
    radius::Float64 
    # fields required for AbstractHamiltonian
    lattice::AbstractLattice
    adjacency::Matrix{Bool}

    StraussHamiltonian(beta::Float64, radius::Float64) = new(beta, radius)
end

function StraussHamiltonian(beta::Float64, radius::Float64, lattice::AbstractLattice)
    sh = StraussHamiltonian(beta, radius)
    lattice!(sh, lattice)
    return sh
end

StraussHamiltonian(beta::Float64, radius::Float64, points::PointSet) = StraussHamiltonian(beta, radius, lattice(points))

mutable struct HardCoreHamiltonian <: AbstractHamiltonian
    radius::Float64
    # fields required for AbstractHamiltonian
    lattice::AbstractLattice
    adjacency::Matrix{Bool}

    HardCoreHamiltonian(radius::Float64) = new(radius)
    
end

function HardCoreHamiltonian(radius::Float64, lattice::AbstractLattice)
    hch = HardCoreHamiltonian(radius)
    lattice!(hch, lattice)
    adjacency_matrix!(hch)
    return hch
end

HardCoreHamiltonian(radius::Float64, points::PointSet) = HardCoreHamiltonian(radius, lattice(points))

function lattice!(sh::StraussHamiltonian, l::AbstractLattice)
    sh.lattice = l
    adjacency_matrix!(sh)
end

function lattice!(hc::HardCoreHamiltonian, l::AbstractLattice)
    hc.lattice = l
    adjacency_matrix!(hc)
end

function local_energy(h::StraussHamiltonian, i::Int)
    return Float64(h.beta*sum(h.adjacency[i, :]))
end

function local_energy(h::HardCoreHamiltonian, i::Int)
    if any(h.adjacency[i, :])
        return Inf
    end
    return 0.0
end

function adjacency_matrix!(h::AbstractHamiltonian)
    points = h.lattice.points
    d = h.lattice.d
    n_points = length(points)
    adjacency = falses(n_points, n_points)

    for i in 1:n_points
        for j in (i+1):n_points
            # if d == 2
            #     dx = points[i][1] - points[j][1]
            #     dy = points[i][2] - points[j][2]
            #     dist_sq = dx*dx + dy*dy
            # elseif d == 3
            #     dx = points[i][1] - points[j][1]
            #     dy = points[i][2] - points[j][2]
            #     dz = points[i][3] - points[j][3]
            #     dist_sq = dx*dx + dy*dy + dz*dz
            # end

            dist_sq = sum((points[i] .- points[j]).^2)
            if dist_sq <= h.radius^2
                adjacency[i, j] = true
                adjacency[j, i] = true
            end
        end
    end

    h.adjacency = adjacency
end


function adjacency_matrix!(h::AbstractHamiltonian, i::Int, new_point::Point)
    n_points = length(h.lattice.points)
    d = h.lattice.d
    points = h.lattice.points
    adjacency = h.adjacency

    # Réinitialiser la ligne et la colonne i
    adjacency[i, :] .= false
    adjacency[:, i] .= false

    # Recalculer les adjacences pour le point i à sa nouvelle position
    for j in 1:n_points
        if j != i
            # if d == 2
            #     dx = new_point[1] - points[j][1]
            #     dy = new_point[2] - points[j][2]
            #     dist_sq = dx*dx + dy*dy
            # elseif d == 3
            #     dx = new_point[1] - points[j][1]
            #     dy = new_point[2] - points[j][2]
            #     dz = new_point[3] - points[j][3]
            #     dist_sq = dx*dx + dy*dy + dz*dz
            # end

            dist_sq = sum((new_point .- points[j]).^2)
            if dist_sq <= h.radius^2
                adjacency[i, j] = true
                adjacency[j, i] = true
            end
        end
    end

end