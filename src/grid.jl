struct Grid <: AbstractLattice
    radius::Int # Domain considered are points on a grid centered and of radius `radius`
    d::Int
    points::Points
    ind::OffsetArray{Int}

    function Grid(radius::Int, d::Int)
        li = OffsetArray{Int}(LinearIndices(tuple(repeat([2radius+1],d)...)), tuple(repeat([-radius:radius],d)...))
        gr = new(radius, d, Points(undef, (2*radius+1)^d), li)
        init!(gr)
        return gr
    end
end

function Base.show(io::IO, g::Grid)
    print(io, "Grid(d=$(g.d),radius=$(g.radius)) with $(length(g)) points")
end

function Base.length(g::Grid)
    return (2*g.radius+1)^g.d
end

function Base.size(g::Grid)
    return tuple(repeat([2*g.radius+1], g.d))
end

dim(g::Grid) = g.d
index(g::Grid) = g.ind

function init!(g::Grid)
    l = length(g)

    if g.d == 2
        for i in 0:(2*g.radius)
            for j in 0:(2*g.radius)
                g.points[(i + 1) + j * (2 * g.radius + 1)] = [Float64(i - g.radius), Float64(j - g.radius)]
            end
        end
    elseif g.d == 3
        for i in 0:(2*g.radius)
            for j in 0:(2*g.radius)
                for k in 0:(2*g.radius)
                    g.points[(i + 1) + j * (2 * g.radius + 1) + k * (2 * g.radius + 1)^2] = [Float64(i - g.radius), Float64(j - g.radius), Float64(k - g.radius)]
                end
            end
        end
    else
        throw(ArgumentError("Dimension d must be 2 or 3"))
    end

end

scale!(g::Grid, scale::Float64) = g.points .*= scale

# SubGrid is just a grid with SubGrid.radius < Grid.radius

