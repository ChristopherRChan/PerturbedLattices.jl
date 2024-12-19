using GeoStats
import Meshes: Vec, Point
using Distributions

g = (10.0, 10.0) # or (10.0, 10.0, 10.0)
d = length(g)
δ = 1.0
q = Normal(0,1)
grid = vertices(CartesianGrid(.-(g), g; dims = Int.(2 .* g ./ δ))) ## vertices() really matters otherwise it is a hexahedron
typeof(grid)
grid isa Vector
vps = PointSet([ v + Vec(rand(q, d)...) for v in grid ])

vps isa PointSet

l2 = sum((vps[1] - vps[2]) .^2)
typeof(l2)
l2.val
Cartesian(1,2)
vps.geoms[1] = Point(1,2)
