struct HamiltonianFamily <: AbstractHamiltonian
    hs::Vector{AbstractHamiltonian}
end

nbparam(h::HamiltonianFamily) = sum([nbparam(h_) for h_=h.hs])
params(h::HamiltonianFamily) = reduce(vcat, [params(h_) for h_=h.hs])
function params!(h::HamiltonianFamily, θ::Vector{Float64})
    inc = 0
    for h_ in h.hs
        if nbparam(h) > 0
            params!(h_, θ[inc + 1:inc+nbparam(h_)])
            inc += nbparam(h_)
        end
    end
end

mutable struct StraussHamiltonian <: AbstractHamiltonian
    β::Float64
    ρ²::Float64 
    # field required for AbstractHamiltonian
    pointset::AbstractPointSet
    StraussHamiltonian(β::Float64, ρ::Float64) = new(β, ρ^2)
end

function StraussHamiltonian(β::Float64, ρ::Float64, ps::AbstractPointSet)
    sh = StraussHamiltonian(β, ρ)
    pointset!(sh, ps)
    return sh
end

StraussHamiltonian(;β::Float64, ρ::Float64) = StraussHamiltonian(β, ρ)

Base.show(io::IO, sh::StraussHamiltonian) = print(io, "Strauss(β=$(sh.β), ρ²=$(sh.ρ²)")

nbparam(sh::StraussHamiltonian) = 1
params(sh::StraussHamiltonian) = sh.β
params!(sh::StraussHamiltonian, β::Float64) = sh.β = β
params!(sh::StraussHamiltonian, θ::Vector{Float64}) = params!(sh, θ[1])

mutable struct HardCoreHamiltonian <: AbstractHamiltonian
    ρ²::Float64 # square of radius
    # fields required for AbstractHamiltonian
    pointset::AbstractPointSet

    HardCoreHamiltonian(ρ::Float64) = new(ρ^2)
    
end

function HardCoreHamiltonian(ρ::Float64, ps::AbstractPointSet)
    hch = HardCoreHamiltonian(ρ)
    pointset!(hch, ps)
    return hch
end

nbparam(hc::HardCoreHamiltonian) = 0
params(::HardCoreHamiltonian) = []

function pointset!(sh::StraussHamiltonian, ps::AbstractPointSet)
    sh.pointset = convert(PointSet, ps)
    update!(sh.pointset)
end

function pointset!(hch::HardCoreHamiltonian, ps::AbstractPointSet)
    hch.pointset = convert(PointSet, ps)
    update!(hch.pointset)
end


S(h::StraussHamiltonian, i::Int) = Σd²(h.pointset, i, h.ρ²)
S(h::StraussHamiltonian, i::Int, point::Point) = Σd²(h.pointset, i, point, h.ρ²)

function localenergy end
hγ = localenergy

localenergy(h::StraussHamiltonian, i::Int) = h.β * S(h,i)
localenergy(h::StraussHamiltonian, i::Int, point::Point) = h.β * S(h,i, point)

localenergy(h::HardCoreHamiltonian, i::Int) = any(d²(h.pointset,i) .<= h.ρ²) ? Inf : 0.0
localenergy(h::HardCoreHamiltonian, i::Int, point::Point) = any(d²(h.pointset, i, point) .<= h.ρ²) ? Inf : 0.0

# function ΔS(h::StraussHamiltonian, i::Int, point::Point)
#     # before move
#     old_S = S(h, i)
#     # after move
#     move!(h.pointset, i, point)
#     return S(h,i) - old_S
# end