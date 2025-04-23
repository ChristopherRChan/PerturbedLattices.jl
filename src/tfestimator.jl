function fit(gpl::GibbsPerturbedLattice, θinit::Float64, νinit::Float64, dataposition::PointSet, dataperturbation::Vector; estimationmethod::Bool=true, optimmethod=NelderMead())
    if estimationmethod
        function  errorfunc(X)
            h = -contrast(gpl, X[1], X[2], dataposition, dataperturbation)
        println(h)
        h
        end
        initial_x = [θinit, νinit]
        res = optimize(errorfunc, initial_x, optimmethod)
        Optim.minimizer(res)
    else
        #TO DO methode variationnel
        (θinit, νinit)
    end
end

function contrast(gpl::GibbsPerturbedLattice , θ::Float64, ν::Float64, dataposition::PointSet, dataperturbation::Vector)
    qpl = 0
    for (i, _) in enumerate(gpl.grid)
        qpl += localenergy([θ], gpl.interactions, i, dataposition[i], dataposition) + ν * g(dataperturbation[i]) + log(quasipartitionfunction(gpl, θ, ν, i, dataposition))
    end
    return(qpl)
end

function quasipartitionfunction(gpl::GibbsPerturbedLattice , θ::Float64, ν::Float64, i::Int64, pts::PointSet; nMC::Int64=10^2)
    d = length(gpl.radius)
    grid = gpl.grid
    partfunc = 0
    for _ in 1:nMC
        h = localenergy([θ], gpl.interactions, i, grid[i]  + Vec(rand(Normal(0, ν), d)...), pts)
        partfunc += exp(-h)/ nMC
    end
    partfunc
end

function localenergy(θ::Vector{Float64}, interactions::Vector{Interaction}, i::Int64, pt::Point, pts::PointSet)
    h=0 
    for (j, inter) in enumerate(interactions)
        h += θ[j] * localenergy(inter, i, pt, pts)
    end
    h
end


function g(perturbation::Vec)
    res = 0.0
    for x in perturbation
        res += x.val ^ 2.0
    end
    res
end