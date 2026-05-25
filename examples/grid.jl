using PerturbedLattices

grid = Grid(10, 2)

dim(grid)
points(grid)

points(grid)[1]

grid[1]
grid[-10,-10]

grid

PerturbedLattices.scale!(grid, .01)
points(grid)
