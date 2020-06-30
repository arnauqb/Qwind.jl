using DrWatson
@quickactivate "Qwind"
using Qwind

using PyPlot
LogNorm = matplotlib.colors.LogNorm
using RegionTrees

# create black hole
M = 1e8 * M_SUN
mdot = 0.5
spin = 0.0
black_hole = BlackHole(M, mdot, spin)

# create radiation
f_uv = 0.85
f_x = 0.15
shielding_density = 2e8
rin = 200.0
radiation = RERadiation(black_hole, f_uv, f_x, shielding_density, rin)

# create simulation grid
r_min = 0.0
r_max = 10000
z_min = 0.0
z_max = 10000
grid = Grid(r_min, r_max, z_min, z_max)

# initial conditions
z0 = 1.0
n0 = 2e8
v0 = 5e7 / C
rfi = 1600
nlines = 40
logspaced = false
initial_conditions = UniformIC(rin, rfi, nlines, z0, n0, v0, logspaced)

# initialize itnegrators
integrators =
    initialize_integrators(radiation, grid, initial_conditions, atol = 1e-8, rtol = 1e-3)

# run integrators
run_integrators!(integrators)

# initialize quadtree
#windtree = initialize_wind_tree(FirstIter())

#run simulation

fig, ax = plot_streamlines(integrators, xh = 2000)
gcf()

# second iteration
windtree = initialize_wind_tree(integrators, n_time = 1000)
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

fig, ax = plot_streamlines(streamlines.lines, yh = 10000, xh = 10000)
gcf()


# third iteration

windtree = initialize_wind_tree(solvers, n_time = 1000)
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

fig, ax = plot_streamlines(streamlines.lines, yh = 10000, xh = 10000)
gcf()

# forth iteration

windtree = initialize_wind_tree(solvers, n_time = 1000)
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

fig, ax = plot_streamlines(streamlines.lines, yh = 10000, xh = 10000)
gcf()
