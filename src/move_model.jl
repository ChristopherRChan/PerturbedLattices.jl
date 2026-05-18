abstract type AbstractMoveModel end

mutable struct GaussianMoveModel <: AbstractMoveModel
    cov::Matrix{Float64}
    d::Int64

    function GaussianMoveModel(cov::Matrix{Float64})
        m = new(cov)
        # verif size(cov,1) == size(cov,2)
        m.d = size(cov,1)
        return m
    end
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


function Base.rand(rng::AbstractRNG, model::GaussianMoveModel)
    return rand(rng, MvNormal(zeros(model.d), model.cov))
end

function Base.rand(rng::AbstractRNG, model::UniformMoveModel)
   return [rand(rng, Uniform(model.bounds[i][1], model.bounds[i][2])) for i in 1:model.d]
end

function Random.rand!(rng::AbstractRNG, ps::PointSet, model::AbstractMoveModel)
    ps.points = [Base.rand(rng, model) for _ in eachindex(ps.lattice.points)]
end

Random.rand!(ps::PointSet, model::AbstractMoveModel) = Random.rand!(Random.default_rng(), ps, model)