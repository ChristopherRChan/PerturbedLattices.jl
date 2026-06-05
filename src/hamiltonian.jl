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
    MultiStraussHamiltonian(β::Vector{Float64}, ρ::Vector{Float64}) = new(β, ρ.^2)
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
    LennardJonesHamiltonian(β₁::Float64, β₂::Float64; d₁::Int=6, d₂::Int=3, R::Float64=1.0) = new(β₁, β₂, d₁, d₂, R)
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
    StraussHamiltonian(β::Float64, ρ::Float64) = new(β, ρ^2)
end

StraussHamiltonian(;β::Float64, ρ::Float64) = StraussHamiltonian(β, ρ)

Base.show(io::IO, sh::StraussHamiltonian) = print(io, "Strauss(h(d)=$(sh.β)× I_[0,$(sqrt(sh.ρ²))](d)")

nbparam(sh::StraussHamiltonian) = 1
params(sh::StraussHamiltonian) = sh.β
params!(sh::StraussHamiltonian, β::Float64) = sh.β = β
params!(sh::StraussHamiltonian, θ::Vector{Float64}) = params!(sh, θ[1])

mutable struct HardCoreHamiltonian <: AbstractHamiltonian
    ρ²::Float64 # square of radius
    HardCoreHamiltonian(ρ::Float64) = new(ρ^2)
    
end

nbparam(hc::HardCoreHamiltonian) = 0
params(::HardCoreHamiltonian) = []


S(h::StraussHamiltonian, i::Int, ps::PointSet) = Σd²(ps, i, h.ρ²)
S(h::StraussHamiltonian, i::Int, point::Point, ps::PointSet) = Σd²(ps, i, point, h.ρ²)

function S(h::MultiStraussHamiltonian, i::Int, point::Point, ps::PointSet)
    s = [Σd²(ps, i, point, ρ²) for ρ²=h.ρ²]
    return s[2:end] .- s[1:end-1]
end
S(h::MultiStraussHamiltonian, i::Int, ps::PointSet) = S(h, i, ps[i], ps)

function S(h::LennardJonesHamiltonian, i::Int, point::Point, ps::PointSet)
    d2 = filter((!=)(0), d²(ps, i, point))
    return [sum((h.R ./ d2 .^ h.d₁)), -sum(h.R ./ d2 .^ h.d₂)]
end
S(h::LennardJonesHamiltonian, i::Int, ps::PointSet) = S(h, i, ps[i], ps)

S(h::ExponentialFamilyHamiltonian, i::Int, point::Point, pl::AbstractPerturbedLatticeModel) = S(h, i, point, pl.pointset)
S(h::ExponentialFamilyHamiltonian, i::Int, pl::AbstractPerturbedLatticeModel) = S(h, i, pl[i], pl)

function localenergy end
hγ = localenergy

localenergy(h::ExponentialFamilyHamiltonian, i::Int, ps::PointSet) = sum(θ(h) .* S(h,i,ps))
localenergy(h::ExponentialFamilyHamiltonian, i::Int, point::Point, ps::PointSet) = sum(S(h,i, point, ps) .* θ(h))
localenergy(h::ExponentialFamilyHamiltonian, i::Int, point::Point, pl::AbstractPerturbedLatticeModel) = localenergy(h, i, point, pl.pointset)
localenergy(h::ExponentialFamilyHamiltonian, i::Int, pl::AbstractPerturbedLatticeModel) = localenergy(h, i, pl[i], pl)


localenergy(h::HardCoreHamiltonian, i::Int, ps::PointSet) = any(d²(ps,i) .<= h.ρ²) ? Inf : 0.0
localenergy(h::HardCoreHamiltonian, i::Int, point::Point, ps::PointSet) = any(d²(ps, i, point) .<= h.ρ²) ? Inf : 0.0