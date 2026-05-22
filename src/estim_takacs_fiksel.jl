## Takacs-Fiksel: TF(f,Γ;θ) = Σₗ DLRₗ(fₗ,Γ;θ)²
# DLRₗ(fₗ,Γ;θ) = Σᵢ (bn(i,xᵢ)fₗ(i,xᵢ,Γᵢᶜ ; θ) - ∫bn(i,y)fₗ(i,y,Γᵢᶜ; θ)Λₙ(i,y,Γᵢᶜ;θ)dy)²
#.           ≈ Σᵢ {(bn(i,Xᵢ)fₗ(i,xᵢ,Γᵢᶜ ; θ) - Σₓ fₗ (i,x,Γᵢᶜ;θ)exp(-θᵀS(i,x,Γᵢᶜ)) / Σₓ exp(-θᵀS(i,xⱼ,Γᵢᶜ)
#  for exponential family:  fₗ(i,xᵢ,Γᵢᶜ ; θ) = Sₗ(i,xᵢ,Γᵢᶜ)
#            = Σₗ Σᵢ {(bn(i,Xᵢ)Sₗ(i,Xᵢ,Γᵢᶜ) - Σₓ Sₗ(i,yⱼ,Γ)exp(-ΣₖθₖSₖ(i,xⱼ,Γᵢᶜ)) / Σₓ exp(-ΣₖθₖSₖ(i,xⱼ,Γᵢᶜ))}


# left term Σf is about point config
# right term ΣΣfΛ is about the inside grid

# default, if hamiltonian or move_model is in TakacsFikselFunctions then it is its default in fns functions
# keyword arg fns is a way to not use the default tf function of hamiltonian ou move_model

# 2 schemes:
# 1. pl is a simulation model
# 2. pl is a model specific to estimation and Γ is its realization from a PointSet supposed to be known (grid structure) eventhough you don't know the link beetwen points and grid points 
# in bith case pl has to be iitialized before TakacsFiksel object

mutable struct TakacsFiksel
    pl::PerturbedLatticeModel
    radius::Int                 # radius defining the inside grid compared with radius of pl.pointset
    fns::Vector{Function}
    # internally
    θ::Vector{Float64}          # final parameters at the end of estimation 
    fcache::Matrix{Float64}     # config cache for left term of DLR
    ΣfΛcache::Matrix{Float64}   # grid cache for right term of DLR
end

function TakacsFiksel(pl::PerturbedLatticeModel, radius::Int; fns::Vector{Function}=Function[])
    θ = zeroes(Float64, nbparam(pl))
    tf = TakacsFiksel(pl, radius, fns, θ, zeros(0, 0), zeros(0, 0))
    init!(tf)
    return tf
end

# Todo: define TakacsFiksel for estimation from pointset
# or define fit(PerturbedLatticeModel, ...) which first define TakacsFiksel object to estimate parameters

function init!(tf::TakacsFiksel) 
    # int fns
    tf.fns = Function[]
    if nbparam(tf.pl.move) > 0
        push!(tf.fns, Base.Fix{1}(tf,tf.pl.move))
    end
    hs = tf.pl.h isa HamiltonianFamily ? tf.pl.h.has : [tf.pl.h]
    for h in hs
       if nbparam(h) > 0
            push!(tf.fns, Base.Fix{1}(tf,h))
        end
    end

    # init cache
end

function prepare_cache!(TakacsFiksel)
    # left cache


    # right cache
end

# one global function to update fields 
function update!(tf::TakacsFiksel; fns::Vector{Function}=Function[], radius::Int=tf.radius)
    if !isempty(fns)
        tf.fns = fns
    end
    if Γ != tf.Γ
        tf.Γ = Γ
    end
    if radius != tf.radius
        tf.radius = radius
    end
end

function Σf_ΣΣfΛ(tf::TakacsFiksel) 
    []
end

function contrast(tf::TakacsFiksel, θ::Vector{Float64})
    params!(tf, θ)
    Σf, ΣΣfΛ = Σf_ΣΣfΛ(tf)
    return sum(((Σf .- ΣΣfΛ).^2))
end

function fit(tf::TakacsFiksel, θ::Vector{Float64})
end

# Takacs Fiksel function test for MoveModel and Hamiltonian
# since Γ is initialized inside plΓ object
function tf_fn(m::GaussianMoveModel, i::Int, x::Point, plΓ::PerturbedLatticeModel)
    sum((x .- points(plΓ)[i]).^2) / 2
end

function tf_fn(h::StraussHamiltonian, i::Int, x::Point)
    
end