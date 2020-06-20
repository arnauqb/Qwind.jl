using Qwind
using Test
using PyPlot
LogNorm = matplotlib.colors.LogNorm

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
radiation = SimpleRadiation(black_hole, f_uv, f_x, shielding_density, r_in)

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
    z_0 = 1e-2,
    log_spaced = true,
)

# initialize quadtree

quadtree = initialize_quadtree(grid)

#run simulation

solvers = []
for line in streamlines.lines
    parameters = Parameters(line, grid, radiation, quadtree)
    solver = initialize_solver(line, parameters)
    push!(solvers, solver)
    run_solver!(solver)
end

function gaussian_smoothing(x, mu, sigma)
    return 1 / (sigma * 2π) * exp(-0.5 * ((x - mu) / sigma)^2)
end
function my_distance(dp, d1, d2)
    ret = d2^2 - ((d1^2 - dp^2 - d2^2) / (2dp))^2
    return sqrt(ret)
end

function initialize_kdtree(solvers)
    solver = solvers[1]
    t_range = 10 .^ range(-3, log10(solver.t), length = 10000)
    dense_solution = solver.sol(t_range)
    r_dense_array = dense_solution[1, :]
    z_dense_array = dense_solution[2, :]
    v_r_dense_array = dense_solution[3, :]
    v_z_dense_array = dense_solution[4, :]
    z_max_array = maximum(z_dense_array) .* ones((length(z_dense_array)))
    widths_array = r_dense_array .* solver.p.line.width_norm
    if length(solvers) >= 2
        for solver in solvers[2:end]
            t_range =
                10 .^ range(-3, log10(solver.sol.t[end-2]), length = 10000)
            dense_solution = solver.sol(t_range)
            r_dense_array = vcat(r_dense_array, dense_solution[1, :])
            z_dense_array = vcat(z_dense_array, dense_solution[2, :])
            v_r_dense_array = vcat(v_r_dense_array, dense_solution[3, :])
            v_z_dense_array = vcat(v_z_dense_array, dense_solution[4, :])
            widths_array = vcat(
                widths_array,
                dense_solution[1, :] .* solver.p.line.width_norm,
            )
            z_max_array = vcat(
                z_max_array,
                maximum(dense_solution[2, :]) .*
                ones((length(dense_solution[2, :]))),
            )
        end
    end
    v_t = @. sqrt(v_r_dense_array^2 + v_z_dense_array^2)
    points = hcat(r_dense_array, z_dense_array)'
    density_dense_array =
        compute_density.(r_dense_array, z_dense_array, v_t, Ref(solver.p.line))
    kdtree = KDTree(points)
    return kdtree,
    density_dense_array,
    widths_array,
    z_max_array,
    r_dense_array,
    z_dense_array
end

function get_density_from_grid(
    kdtree,
    r,
    z,
    densities,
    widths,
    r_dense_array,
    z_dense_array,
)
    idcs, distances = knn(kdtree, [r, z], 2)
    p1 = [r_dense_array[idcs[1]], z_dense_array[idcs[1]]]
    d1 = distances[1]
    p2 = [r_dense_array[idcs[2]], z_dense_array[idcs[2]]]
    d2 = distances[2]
    dp = evaluate(Euclidean(), p1, p2)
    distance = my_distance(dp, d1, d2)
    width = sum(widths[idcs]) / 4
    density = sum(densities[idcs]) / 2
    if distance <= width
        ret = density #max(density, density_grid[i,j])
    else
        ret = max(density * exp(-(distance - width)^2 / 2), 1e2)
    end
    return ret
end

kdtree, densities, widths, z_array, r_dense_array, z_dense_array =
    initialize_kdtree(solvers)
MIN_R = 0
MIN_Z = 0
MAX_R = 1000
MAX_Z = 1000#600#0.21
r_range = range(MIN_R, MAX_R, length = convert(Int, 500))
z_range = range(MIN_Z, MAX_Z, length = convert(Int, 501))
density_grid = 1e2 .* ones((length(r_range), length(z_range)))
for (i, r) in enumerate(r_range)
    for (j, z) in enumerate(z_range)
        density_grid[i,j] = get_density_from_grid(kdtree, r, z, densities, widths, r_dense_array, z_dense_array)
    end
end
fig, ax = plt.subplots()
cm = ax.pcolormesh(r_range, z_range, density_grid', norm = LogNorm())
plt.colorbar(cm, ax = ax)
#ax.plot(r_dense_array, z_dense_array, "o", markersize=0.1)
plot_streamlines(streamlines.lines, fig, ax)
#for i in 1:length(streamlines)
#    idx = (i-1)*1000+1:i*1000
#    ax.plot(r_dense_array[idx], z_dense_array[idx], "o", markersize=1)
#end
ax.set_xlim(MIN_R, MAX_R)
ax.set_ylim(MIN_Z, MAX_Z)
#ax.set_yscale("log")
gcf()
#xlims!(0,1000)
#ylims!(0,100)
