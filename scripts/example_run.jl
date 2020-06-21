using Qwind
using Test
using PyPlot
LogNorm = matplotlib.colors.LogNorm
using RegionTrees
#pyplot()

# create black hole
M = 1e8 * M_SUN
mdot = 0.5
spin = 0.0
black_hole = BlackHole(M, mdot, spin)

# create radiation
f_uv = 0.85
f_x = 0.15
shielding_density = 2e8
r_in = 20.0
radiation = RERadiation(black_hole, f_uv, f_x, shielding_density, r_in)

# create simulation grid
r_min = 0.0
r_max = 10000
z_min = 0.0
z_max = 10000
vacuum_density = 1e2
grid = Grid(r_min, r_max, z_min, z_max)

# initialize streamlines
r_fi = 1600
n_lines = 80
velocity_generator(r_0) = 5e7
density_generator(r_0) = 2e8
streamlines = Streamlines(
    r_in,
    r_fi,
    n_lines,
    velocity_generator,
    density_generator,
    black_hole.M;
    z_0 = 1.0,
    log_spaced = false,
)

# initialize quadtree
windtree = initialize_wind_tree(FirstIter())

#run simulation

# first iteration
solvers = []
for line in streamlines.lines
    parameters = Parameters(line, grid, radiation, windtree)
    solver = initialize_solver(line, parameters)
    run_solver!(solver)
    push!(solvers, solver)
end

fig, ax = plot_streamlines(streamlines.lines, xh=2000)
gcf()

# second iteration
windtree = initialize_wind_tree(solvers, n_time=1000)
refine_quadtree!(windtree)

r_range = range(0, 1000, length = 500)
z_range = range(0, 1000, length = 500)

density_grid = zeros((length(r_range), length(z_range)))

for (i, r) in enumerate(r_range)
    for (j, z) in enumerate(z_range)
        density_grid[i, j] = findleaf(windtree.quadtree, [r, z]).data
    end
end

fig, ax = plt.subplots()
cm = ax.pcolormesh(r_range, z_range, density_grid', norm = LogNorm())
plt.colorbar(cm, ax = ax)
gcf()

streamlines = Streamlines(
    r_in,
    r_fi,
    n_lines,
    velocity_generator,
    density_generator,
    black_hole.M;
    z_0 = 1.0,
    log_spaced = false,
)

solvers = []
for line in streamlines.lines
    parameters = Parameters(line, grid, radiation, windtree)
    solver = initialize_solver(line, parameters)
    run_solver!(solver)
    push!(solvers, solver)
end

fig, ax = plot_streamlines(streamlines.lines, yh=10000, xh=10000)
gcf()


# third iteration

windtree = initialize_wind_tree(solvers, n_time=1000)
refine_quadtree!(windtree)

r_range = range(0, 1000, length = 500)
z_range = range(0, 1000, length = 500)

density_grid = zeros((length(r_range), length(z_range)))

for (i, r) in enumerate(r_range)
    for (j, z) in enumerate(z_range)
        density_grid[i, j] = findleaf(windtree.quadtree, [r, z]).data
    end
end

fig, ax = plt.subplots()
cm = ax.pcolormesh(r_range, z_range, density_grid', norm = LogNorm())
plt.colorbar(cm, ax = ax)
gcf()

streamlines = Streamlines(
    r_in,
    r_fi,
    n_lines,
    velocity_generator,
    density_generator,
    black_hole.M;
    z_0 = 1.0,
    log_spaced = false,
)

solvers = []
for line in streamlines.lines
    parameters = Parameters(line, grid, radiation, windtree)
    solver = initialize_solver(line, parameters)
    run_solver!(solver)
    push!(solvers, solver)
end

fig, ax = plot_streamlines(streamlines.lines, yh=10000, xh=10000)
gcf()

# forth iteration

windtree = initialize_wind_tree(solvers, n_time=1000)
refine_quadtree!(windtree)

r_range = range(0, 1000, length = 500)
z_range = range(0, 1000, length = 500)

density_grid = zeros((length(r_range), length(z_range)))

for (i, r) in enumerate(r_range)
    for (j, z) in enumerate(z_range)
        density_grid[i, j] = findleaf(windtree.quadtree, [r, z]).data
    end
end

fig, ax = plt.subplots()
cm = ax.pcolormesh(r_range, z_range, density_grid', norm = LogNorm())
plt.colorbar(cm, ax = ax)
gcf()

streamlines = Streamlines(
    r_in,
    r_fi,
    n_lines,
    velocity_generator,
    density_generator,
    black_hole.M;
    z_0 = 1.0,
    log_spaced = false,
)

solvers = []
for line in streamlines.lines
    parameters = Parameters(line, grid, radiation, windtree)
    solver = initialize_solver(line, parameters)
    run_solver!(solver)
    push!(solvers, solver)
end

fig, ax = plot_streamlines(streamlines.lines, yh=10000, xh=10000)
gcf()
