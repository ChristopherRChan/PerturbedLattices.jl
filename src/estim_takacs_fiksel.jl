## Takacs-Fiksel: TF(f,Γ;θ) = Σₗ DLRₗ(fₗ,Γ;θ)²
# DLRₗ(fₗ,Γ;θ) = Σᵢ (bn(i,xᵢ)fₗ(i,xᵢ,Γᵢᶜ ; θ) - ∫bn(i,y)fₗ(i,y,Γᵢᶜ; θ)Λₙ(i,y,Γᵢᶜ;θ)dy)²
#            ≈ Σᵢ {(bn(i,Xᵢ)fₗ(i,xᵢ,Γᵢᶜ ; θ) - Σₓ fₗ (i,x,Γᵢᶜ;θ)exp(-θᵀS(i,x,Γᵢᶜ)) / Σₓ exp(-θᵀS(i,xⱼ,Γᵢᶜ)
# for exponential family:  fₗ(i,xᵢ,Γᵢᶜ ; θ) = Sₗ(i,xᵢ,Γᵢᶜ)
#            = Σₗ Σᵢ {(bn(i,Xᵢ)Sₗ(i,Xᵢ,Γᵢᶜ) - Σₓ Sₗ(i,yⱼ,Γ)exp(-ΣₖθₖSₖ(i,xⱼ,Γᵢᶜ)) / Σₓ exp(-ΣₖθₖSₖ(i,xⱼ,Γᵢᶜ))}


# left term Σf is about point config
# right term ΣΣfΛ is about the inside grid

# default, if hamiltonian or move_model is in TakacsFikselFunctions then it is its default in f functions
# keyword arg f is a way to not use the default tf function of hamiltonian ou move_model

# 2 schemes:
# 1. pl is a simulation model
# 2. pl is a model specific to estimation and Γ is its realization from a PointSet supposed to be known (grid structure) eventhough you don't know the link beetwen points and grid points 
# in bith case pl has to be iitialized before TakacsFiksel object

mutable struct TakacsFiksel
    pl::PerturbedLatticeModel
    radius::Int                 # radius defining the inside grid compared with radius of pl.pointset
    f::Vector{Function}
    # internally
    θ::Vector{Float64}          # final parameters at the end of estimation 
    fcache::Matrix{Float64}     # config cache for left term of DLR
    ΣfΛcache::Array{Float64}   # grid cache for right term of DLR
    gridQ::Grid
end

function TakacsFiksel(pl::PerturbedLatticeModel, radius::Int; f::Vector{Function}=Function[], gridQ::Grid=Grid(2,2))
    θ = zeros(Float64, nbparam(pl))
    tf = TakacsFiksel(pl, radius, f, θ, zeros(0, 0), zeros(0, 0), gridQ)
    init!(tf)
    return tf
end

# Todo: define TakacsFiksel for estimation from pointset
# or define fit(PerturbedLatticeModel, ...) which first define TakacsFiksel object to estimate parameters

function init!(tf::TakacsFiksel) 
    # int f
    tf.f = Function[]
    if nbparam(tf.pl.move) > 0
        push!(tf.f, Base.Fix{1}(fₗ,tf.pl.move))
    end
    hs = tf.pl.h isa HamiltonianFamily ? tf.pl.h.hs : [tf.pl.h]
    for h in hs
       if nbparam(h) > 0
            push!(tf.f, Base.Fix{1}(fₗ,h))
        end
    end
end

function cache!(tf::TakacsFiksel; nQ::Int=50, ρQ::Float64=3.0)
    d = tf.pl.pointset.lattice.d
    subind = lattice(tf.pl).ind[repeat([-tf.radius:tf.radius],d)...]
    
    # left cache
    tf.fcache = Matrix{Float64}(undef, length(subind), length(tf.f))
    pts = points(tf.pl)
    for (i, si) in enumerate(subind) # si is the index in the subgrid
        for l in eachindex(length(tf.f))
            tf.fcache[i, l] = tf.f[l](si, pts[si], tf.pl)
        end
    end 

    # right cache
    tf.gridQ = Grid(nQ , d)
    δ = ρQ / nQ
    scale!(tf.gridQ, δ)
    tf.ΣfΛcache = Array{Float64}(undef, length(subind), length(tf.gridQ) ,length(tf.f))
    lpts = points(lattice(tf.pl))
    for (i, si) in enumerate(subind)
        for k in eachindex(points(tf.gridQ))
            pt = lpts[si] .+ tf.gridQ[k]
            for l in eachindex(length(tf.f))
                tf.ΣfΛcache[i, k, l] = tf.f[l](si, pt, tf.pl)
            end 
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
fₗ(m::GaussianMoveModel, i::Int, x::Point, plΓ::PerturbedLatticeModel) = sum((x .- points(lattice(plΓ))[i]).^2) / 2

fₗ(h::StraussHamiltonian, i::Int, x::Point, plΓ::PerturbedLatticeModel) = S(h, i, point)