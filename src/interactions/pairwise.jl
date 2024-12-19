struct PairwiseInteraction <: Interaction
    φ::Function # on l²
end

pairwise(φ::Function) = PairwiseInteraction(φ ∘ √)

function localenergy(pair::PairwiseInteraction, i::Int64, pt::Point, points::PointSet)
    locen = 0.0
    for j in eachindex(points)
        if j != i
            l² = sum((pt - points[j]).^2).val
            locen += pair.φ(l²)
        end
    end
    return locen
end

function moveenergy(pair::PairwiseInteraction, i::Int64, newpoint::Point, points::PointSet)
     return localenergy(pair, i, newpoint, points) - localenergy(pair, i, points[i], points)
end