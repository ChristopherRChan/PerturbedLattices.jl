## Takacs-Fiksel: TF(f,Γ;θ) = Σₗ DLRₗ(fₗ,Γ;θ)²
# DLRₗ(fₗ,Γ;θ) = Σᵢ (bn(i,xᵢ)fₗ(i,xᵢ,Γᵢᶜ ; θ) - ∫bn(i,y)fₗ(i,y,Γᵢᶜ; θ)Λₙ(i,y,Γᵢᶜ;θ)dy)²
#            ≈ Σᵢ {(bn(i,Xᵢ)fₗ(i,xᵢ,Γᵢᶜ ; θ) - Σₓ fₗ (i,x,Γᵢᶜ;θ)exp(-θᵀS(i,x,Γᵢᶜ)) / Σₓ exp(-θᵀS(i,xⱼ,Γᵢᶜ)
# for exponential family:  fₗ(i,xᵢ,Γᵢᶜ ; θ) = Sₗ(i,xᵢ,Γᵢᶜ)
#            = Σₗ Σᵢ {(bn(i,Xᵢ)Sₗ(i,Xᵢ,Γᵢᶜ) - Σₓ Sₗ(i,yⱼ,Γᵢᶜ)exp(-ΣₖθₖSₖ(i,xⱼ,Γᵢᶜ)) / Σₓ exp(-ΣₖθₖSₖ(i,xⱼ,Γᵢᶜ))}


# left term Σf is about point config
# right term ΣΣfΛ is about the inside grid

# default, if hamiltonian or move_model is in TakacsFikselFunctions then it is its default in f functions
# keyword arg f is a way to not use the default tf function of hamiltonian ou move_model

# 2 schemes:
# 1. pl is a simulation model
# 2. pl is a model specific to estimation and Γ is its realization from a PointSet supposed to be known (grid structure) eventhough you don't know the link beetwen points and grid points 
# in bith case pl has to be iitialized before TakacsFiksel object

mutable struct TakacsFiksel <: EstimationMethod
    pl::PerturbedLatticeModel
    radius::Int                 # radius defining the inside grid compared with radius of pl.pointset
    f::Vector{Function}         # f[l] returns Vector{Float64} corresponding to [fₗ(...) for l=1:nθ(pl.h)]
    # internally
    θ::Vector{Float64}          # final parameters at the end of estimation 
    fcache::Matrix{Float64}     # config cache for left term of DLR
    ΣfΛcache::Array{Float64}   # grid cache for right term of DLR
    gridQ::Grid
    subind::Array{Int}
end

function TakacsFiksel(pl::PerturbedLatticeModel, radius::Int; f::Vector{Function}=Function[], nQ::Int=5, ρQ::Float64=1.5)
    θ = zeros(Float64, nbparam(pl))
    # grid for approx integral on Q
    d = dim(pl)
    gridQ = Grid(nQ , d)
    δ = ρQ / nQ
    scale!(gridQ, δ)
    subind = lattice(pl).ind[repeat([-radius:radius],d)...]
    tf = TakacsFiksel(pl, radius, f, θ, zeros(0, 0), zeros(0, 0), gridQ, subind)
    init!(tf)
    return tf
end

Base.show(io::IO, tf::TakacsFiksel) = print(io, "TakacsFiksel(Γ=$(tf.pl))")

nbparam(tf::TakacsFiksel) = nbparam(tf.pl)
params(tf::TakacsFiksel) = params(tf.pl)
params!(tf::TakacsFiksel, θ::Vector{Float64}) = params!(tf.pl, θ)

points(tf::TakacsFiksel) = points(tf.pl)

# Todo: define TakacsFiksel for estimation from pointset
# or define fit(PerturbedLatticeModel, ...) which first define TakacsFiksel object to estimate parameters

function init!(tf::TakacsFiksel) 
    # int f
    tf.f = Function[]
    if nbparam(tf.pl.move) > 0
        push!(tf.f, Base.Fix1(fₗ,tf.pl.move))
    end
    hs = tf.pl.h isa HamiltonianFamily ? tf.pl.h.hs : [tf.pl.h]
    for h in hs
       if nbparam(h) > 0
            push!(tf.f, Base.Fix1(fₗ,h))
        end
    end
end

function cache!(tf::TakacsFiksel)
    # left cache
    tf.fcache = Matrix{Float64}(undef, length(tf.subind), nθ(tf))
    pts = points(tf.pl)
    for (i, si) in enumerate(tf.subind) # si is the index in the subgrid
        tf.fcache[i, :] = vcat([tf.f[l](si, pts[si], tf.pl) for l in eachindex(tf.f)]...)
    end 

    # right cache
    
    tf.ΣfΛcache = Array{Float64}(undef, length(tf.subind), length(tf.gridQ) ,nθ(tf))
    lpts = points(lattice(tf.pl))
    for (i, si) in enumerate(tf.subind)
        for k in eachindex(points(tf.gridQ))
            pt = lpts[si] .+ tf.gridQ[k]
            tf.ΣfΛcache[i, k, :] = vcat([tf.f[l](si, pt, tf.pl) for l in eachindex(tf.f)]...)
        end
    end
end

# one global function to update fields 
function update!(tf::TakacsFiksel; f::Vector{Function}=Function[], radius::Int=tf.radius)
    if !isempty(f)
        tf.f = f
    end
    if Γ != tf.Γ
        tf.Γ = Γ
    end
    if radius != tf.radius
        tf.radius = radius
    end
end

contrast(tf::TakacsFiksel) = sum(
    (
        transpose(fill(1.0, length(tf.subind))) * 
            (
                tf.fcache -  
                    [ sum((exp.(-tf.ΣfΛcache[i, : , :]  * θ(tf)) .* tf.ΣfΛcache[i, : , l])) / sum(exp.(-tf.ΣfΛcache[i, : , :]  * θ(tf))) for i=eachindex(tf.subind), l=eachindex(θ(tf))]
            ) / length(tf.subind)
    ) .^2
)

function contrast(tf::TakacsFiksel,::Type{Vector})
    res = zeros(length(tf.f))
    for i=eachindex(tf.subind)
        for l=eachindex(θ(tf))
            res[l] += tf.fcache[i, l] 
            ΣS, Σ = 0.0, 0.0
            for k=eachindex(points(tf.gridQ))
                tmp = sum(exp.(-tf.ΣfΛcache[i, k , :]'  * θ(tf)))
                Σ += tmp
                ΣS += tmp * tf.ΣfΛcache[i, k , l]
            end
            res[l] -= ΣS / Σ
        end
    end
    return sum((res ./ length(tf.subind)) .^2)
end

function contrast(tf::TakacsFiksel, θ::Vector{Float64})
    θ!(tf, θ)
    return contrast(tf) #, Vector)
end

function fit!(tf::TakacsFiksel, θ₀::Vector{Float64} = θ(tf))
    cache!(tf)
    θ!(tf, θ₀) 
    result = optimize(Base.Fix1(contrast, tf), θ₀, NelderMead())
    tf.θ = Optim.minimizer(result)
    result
end

# Takacs Fiksel function test for Move and Hamiltonian
# since Γ is initialized inside plΓ object
fₗ(m::GaussianMove, i::Int, x::Point, plΓ::PerturbedLatticeModel) = sum((x .- points(lattice(plΓ))[i]).^2) / 2

fₗ(h::ExponentialFamilyHamiltonian, i::Int, x::Point, plΓ::PerturbedLatticeModel) = S(h, i, x, plΓ)

