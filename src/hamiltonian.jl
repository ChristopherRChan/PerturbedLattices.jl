function pointset!(h::AbstractHamiltonian, ps::AbstractPointSet)
    h.pointset = convert(PointSet, ps)
    update!(h.pointset)
end

struct HamiltonianFamily <: ExponentialFamilyHamiltonian
    hs::Vector{ExponentialFamilyHamiltonian}
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

mutable struct MultiStraussHamiltonian <: ExponentialFamilyHamiltonian
    β::Vector{Float64}
    ρ²::Vector{Float64} # ρ²[1]=HC and ρ²[end]=Range (length(ρ²) = length(β) + 1)
    # field required for AbstractHamiltonian
    pointset::AbstractPointSet
    MultiStraussHamiltonian(β::Vector{Float64}, ρ::Vector{Float64}) = new(β, ρ.^2)
end

function MultiStraussHamiltonian(β::Vector{Float64}, ρ::Vector{Float64}, ps::AbstractPointSet)
    h = MultiStraussHamiltonian(β, ρ)
    pointset!(h, ps)
    return h
end

MultiStraussHamiltonian(;β::Vector{Float64}, ρ::Vector{Float64}) = MultiStraussHamiltonian(β, ρ)

function Base.show(io::IO, h::MultiStraussHamiltonian)
    print(io, "MultiStrauss(h(d)=")
    if h.ρ²[1] > 0
        print(io, "(∞) × I_[0,$(sqrt(h.ρ²[1]))](d)")
        if length(h.ρ²) > 1
            print(io, " + ")
        end
    end
    if length(h.ρ²) > 1
        print(io, join(["$(h.β[i])× I_[$(sqrt(h.ρ²[i])),$(sqrt(h.ρ²[i+1]))](d)" for i=eachindex(h.β)], " * "))
    end
    print(io, ")")
end

nbparam(h::MultiStraussHamiltonian) = length(h.β)
params(h::MultiStraussHamiltonian) = h.β
params!(h::MultiStraussHamiltonian, β::Vector{Float64}) = h.β = β



mutable struct LennardJonesHamiltonian <: ExponentialFamilyHamiltonian
    β₁::Float64
    β₂::Float64
    d₁::Int
    d₂::Int
    R::Float64
    # field required for AbstractHamiltonian
    pointset::AbstractPointSet
    LennardJonesHamiltonian(β₁::Float64, β₂::Float64; d₁::Int=6, d₂::Int=3, R::Float64=1.0) = new(β₁, β₂, d₁, d₂, R)
end

function LennardJonesHamiltonian(β₁::Float64, β₂::Float64, ps::AbstractPointSet; d₁::Int=12, d₂::Int=6, R::Float64=1.0)
    h = LennardJonesHamiltonian(β₁, β₂, d₁=d₁, d₂=d₂, R=R)
    pointset!(h, ps)
    return h
end

function Base.show(io::IO, h::LennardJonesHamiltonian)
    print(io, "LennardJones(h(d)=$(h.β₁)×($(h.R)/d²)^$(h.d₁)) - $(h.β₂)×($(h.R)/d²)^$(h.d₂)))")
end

nbparam(h::LennardJonesHamiltonian) = 2
params(h::LennardJonesHamiltonian) = [h.β₁, h.β₂]
params!(h::LennardJonesHamiltonian, β::Vector{Float64}) = (h.β₁, h.β₂) = β


mutable struct StraussHamiltonian <: ExponentialFamilyHamiltonian
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

Base.show(io::IO, sh::StraussHamiltonian) = print(io, "Strauss(h(d)=$(sh.β)× I_[0,$(sqrt(sh.ρ²))](d)")

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


S(h::StraussHamiltonian, i::Int) = Σd²(h.pointset, i, h.ρ²)
S(h::StraussHamiltonian, i::Int, point::Point) = Σd²(h.pointset, i, point, h.ρ²)

function S(h::MultiStraussHamiltonian, i::Int, point::Point)
    s = [Σd²(h.pointset, i, point, ρ²) for ρ²=h.ρ²]
    return s[2:end] .- s[1:end-1]
end
S(h::MultiStraussHamiltonian, i::Int) = S(h, i, h.pointset[i])

function S(h::LennardJonesHamiltonian, i::Int, point::Point)
    [sum((h.R ./ (filter((!=)(0), d²(h.pointset, i, point))) .^ h.d₁)), -sum(h.R ./filter((!=)(0), d²(h.pointset, i, point)) .^ h.d₂)]
end
S(h::LennardJonesHamiltonian, i::Int) = S(h, i, h.pointset[i])


function localenergy end
hγ = localenergy

localenergy(h::ExponentialFamilyHamiltonian, i::Int) = sum(θ(h) .* S(h,i))
localenergy(h::ExponentialFamilyHamiltonian, i::Int, point::Point) = sum(S(h,i, point) .* θ(h)) 

localenergy(h::HardCoreHamiltonian, i::Int) = any(d²(h.pointset,i) .<= h.ρ²) ? Inf : 0.0
localenergy(h::HardCoreHamiltonian, i::Int, point::Point) = any(d²(h.pointset, i, point) .<= h.ρ²) ? Inf : 0.0

# function ΔS(h::StraussHamiltonian, i::Int, point::Point)
#     # before move
#     old_S = S(h, i)
#     # after move
#     move!(h.pointset, i, point)
#     return S(h,i) - old_S
# end