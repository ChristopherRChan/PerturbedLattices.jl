mutable struct GaussianMove <: AbstractMove
    Σ::Matrix{Float64}
    d::Int64

    function GaussianMove(cov::Matrix{Float64})
        m = new(cov)
        # verif size(cov,1) == size(cov,2)
        m.d = size(cov,1)
        return m
    end
end

GaussianMove(σ²::Float64, d::Int) = GaussianMove(diagm(repeat([σ²],d)))
GaussianMove(;σ²::Float64=1.0, d::Int=2) = GaussianMove(σ², d)

function Base.show(io::IO, m::GaussianMove)
    print(io, "GaussianMove(σ²=$(θ(m)[1]), d=$(m.d))")
end

# param is 1/σ² = σ⁻²
nbparam(m::GaussianMove) = 1
params(m::GaussianMove) = [1 / m.Σ[1,1]]
params!(m::GaussianMove, σ⁻²::Float64) = m.Σ[diagind(m.Σ)] .= σ⁻²
params!(m::GaussianMove, θ::Vector{Float64}) = params!(m,1/θ[1])


function Base.rand(rng::AbstractRNG, model::GaussianMove)
    return rand(rng, MvNormal(zeros(model.d), model.Σ))
end

mutable struct UniformMove <: AbstractMove
    bounds::Vector{Vector{Float64}}
    d::Int64

    function UniformMove(bounds::Vector{Vector{Float64}})
        m = new(bounds)
        m.d = length(bounds)
        return m
    end
end

function Base.show(io::IO, m::UniformMove)
    print(io, "UniformMove([$(m.bounds[1][1]), $(m.bounds[1][2])]×[$(m.bounds[2][1]), $(m.bounds[2][2])]")
end

nbparam(m::UniformMove) = 0
params(::UniformMove) = []

function Base.rand(rng::AbstractRNG, model::UniformMove)
   return [rand(rng, Uniform(model.bounds[i][1], model.bounds[i][2])) for i in 1:model.d]
end

function Random.rand!(rng::AbstractRNG, ps::PointSet, model::AbstractMove)
    ps.points = [Base.rand(rng, model) for _ in eachindex(ps.lattice.points)]
end

# rand for PointSet

Random.rand!(ps::PointSet, model::AbstractMove) = Random.rand!(Random.default_rng(), ps, model)