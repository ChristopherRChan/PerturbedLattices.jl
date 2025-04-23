function GeoStatsProcesses.randsingle(rng::AbstractRNG, gpl::GibbsPerturbedLattice, g) #g=radius of box
    cpt = 0
    d = length(gpl.radius)
    n = length(gpl.points)
    for _ in 1:gpl.nsim
        i = rand(rng, 1:n)
        newpt = gpl.grid[i] + Vec(rand(rng, gpl.q, d)...)
        if rand() <= exp( -moveenergy(gpl, i, newpt) )
            cpt += 1
            gpl.points.geoms[i] = newpt
        end
    end
    println(cpt, " moves over ", gpl.nsim)
    return gpl.points
end

function moveenergy(gpl::GibbsPerturbedLattice, i::Int64, pt::Point)
    en = 0.0
    for interaction in gpl.interactions
        en += moveenergy(interaction, i, pt, gpl.points)
    end
    return en
end