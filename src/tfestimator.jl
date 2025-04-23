function fit(gpl::GibbsPerturbedLattice, θinit::Float64, νinit::Float64, dataposition::PointSet, dataperturbation::Vector; estimationmethod::Bool=true, optimmethod=NelderMead())
    if estimationmethod
        function  errorfunc(X)
            -contrast(gpl, X[1], X[2], dataposition, dataperturbation)
        end
        initial_x =[θinit, νinit]
        res = optimize(errorfunc, initial_x, optimmethod)
        Optim.minimizer(res)
    else
        #TO DO methode variationnel
        (θinit, νinit)
    end
end

function contrast(gpl::GibbsPerturbedLattice , θ::Float64, ν::Float64, dataposition::PointSet, dataperturbation::Vector; nMC::Int64=10^4)
    qpl = 0
    for i in 1:length(gpl.grid)
        qpl = + localenergy([θ], gpl.interactions, i, dataposition[i], dataposition) + ν*g(dataperturbation[i]) + log(quasipartitionfunction(gpl, θ, ν, i, dataposition))
    end
    return(qpl)
end

function quasipartitionfunction(gpl::GibbsPerturbedLattice , θ::Float64, ν::Float64, i::Int64, pts::PointSet; nMC::Int64=10^4)
    d = length(gpl.radius)
    grid = gpl.grid
    partfunc = 0
        for j in 1:nMC
            h= localenergy([θ], gpl.interactions, i, grid[i]  + Vec(rand.(Normal(0, ν), d)), pts)
            partfunc = +exp(-h)/ nMC
        end
end

function localenergy(θ::Vector{Float64}, interactions::Vector{Interaction}, i::Int64, pt::Point, pts::PointSet)
    h=0 
    for j in 1:length(interactions)
        h =+ θ[j]*localenergy(interactions[j], i, pt, pts)
    end
    return(h)
end


function g(perturbation::Meshes.Vec)
    res =0
    for x in perturbation
        res =+ x^2
    end
    return(res)
end