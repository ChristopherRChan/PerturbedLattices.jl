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


S(h::StraussHamiltonian, i::Int) = sum(square_distance(h.pointset)[i, :] .<= h.radius^2)
local_energy(h::StraussHamiltonian, i::Int) = h.beta * S(h,i)
local_energy(h::HardCoreHamiltonian, i::Int) = any(h.adjacency[i, :] .<= h.radius^2) ? Inf : 0.0

function ΔS(h::StraussHamiltonian, i::Int, point::Point)
    # before move
    old_S = S(h, i)
    # after move
    move!(h.pointset, i, point)
    return S(h,i) - old_S
end