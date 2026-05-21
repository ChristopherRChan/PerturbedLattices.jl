abstract type AbstractMoveModel end

mutable struct GaussianMoveModel <: AbstractMoveModel
    Σ::Matrix{Float64}
    d::Int64

    function GaussianMoveModel(cov::Matrix{Float64})
        m = new(cov)
        # verif size(cov,1) == size(cov,2)
        m.d = size(cov,1)
        return m
    end
end

GaussianMoveModel(σ²::Float64, d::Int) = GaussianMoveModel(diagm(repeat([σ²],d)))
nbparam(m::GaussianMoveModel) = 1
params(m::GaussianMoveModel) = [m.Σ[1,1]]
params!(m::GaussianMoveModel, σ²::Float64) = m.Σ[diagind(m.Σ)] .= σ²
params!(m::GaussianMoveModel, θ::Vector{Float64}) = params!(m,θ[1])


function Base.rand(rng::AbstractRNG, model::GaussianMoveModel)
    return rand(rng, MvNormal(zeros(model.d), model.Σ))
end

mutable struct UniformMoveModel <: AbstractMoveModel
    bounds::Vector{Vector{Float64}}
    d::Int64

    function UniformMoveModel(bounds::Vector{Vector{Float64}})
        m = new(bounds)
        m.d = length(bounds)
        return m
    end
end

nbparam(m::UniformMoveModel) = 0
params(::UniformMoveModel) = []

function Base.rand(rng::AbstractRNG, model::UniformMoveModel)
   return [rand(rng, Uniform(model.bounds[i][1], model.bounds[i][2])) for i in 1:model.d]
end

function Random.rand!(rng::AbstractRNG, ps::PointSet, model::AbstractMoveModel)
    ps.points = [Base.rand(rng, model) for _ in eachindex(ps.lattice.points)]
end

# rand for PointSet

Random.rand!(ps::PointSet, model::AbstractMoveModel) = Random.rand!(Random.default_rng(), ps, model)