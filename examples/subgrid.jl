using PerturbedLattices

g = Grid(20, 2)

g[-10,10]

g.ind[-10:10,-10:10]

vi = @view g.ind[-10:10,-10:10]

vi[1]
vi[end]