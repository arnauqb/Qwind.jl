using DrWatson
@quickactivate "Qwind"
using YAML
using PyCall, RegionTrees
using PyPlot
LogNorm = matplotlib.colors.LogNorm
@pyimport scipy.interpolate as interp
using Qwind

function plot_density_grid(
    quadtree::Cell;
    xl = 0,
    xh = 5000,
    yl = 0,
    yh = 5000,
    nr = 500,
    nz = 501,
    vmin = nothing,
    vmax = nothing,
    fig = nothing,
    ax = nothing,
)
    if fig === nothing || ax === nothing
        fig, ax = plt.subplots()
    end
    grid = zeros(Float64, nr, nz)
    r_range = range(xl, stop = xh, length = nr)
    z_range = range(yl, stop = yh, length = nz)
    for (i, r) in enumerate(r_range[1:(end - 1)])
        for (j, z) in enumerate(z_range[1:(end - 1)])
            grid[i, j] = get_density(quadtree, r, z)
        end
    end
    cm = ax.pcolormesh(r_range, z_range, grid', norm = LogNorm())
    plt.colorbar(cm, ax = ax)
    ax.set_xlabel("r [Rg]")
    ax.set_ylabel("z [Rg]")
    ax.set_xlim(xl, xh)
    ax.set_ylim(yl, yh)
    return fig, ax
end

config = YAML.load_file("scripts/config_uniform_ic.yaml");

black_hole = BlackHole(config);

radiation = @eval $(Symbol(config["radiation"]["mode"]))(black_hole, config);

radiative_transfer =
    @eval $(Symbol(config["radiative_transfer"]["mode"]))(radiation, config);

grid = Rectangular(config);

initial_conditions =
    @eval $(Symbol(config["initial_conditions"]["mode"]))(radiation, black_hole, config);

integrators = initialize_integrators(
    radiative_transfer,
    grid,
    initial_conditions,
    atol = config["integrator"]["atol"],
    rtol = config["integrator"]["rtol"],
);

run_integrators!(integrators);

#r, z, zmax, z0, line_width, density = get_dense_solution_from_integrators(integrators, 1000);

kdtree = create_wind_kdtree(integrators, 100000);

smoothed_grid = SmoothedGrid(kdtree, nr = 500, nz = 501, kernel_size = 1);

fig, ax = plt.subplots()
cm = ax.pcolormesh(
    smoothed_grid.r_range,
    smoothed_grid.z_range,
    smoothed_grid.grid',
    norm = LogNorm(1e2, 1e9),
)
plot_streamlines(integrators, fig, ax)
plt.colorbar(cm, ax = ax)
ax.set_xlim(0, 100)
ax.set_ylim(0, 100)


quadtree = create_and_refine_quadtree(
    smoothed_grid,
    Rg = black_hole.Rg,
    atol = 1e-3,
    rtol = 1e-3,
    cell_min_size = 1e-6,
)


fig, ax = plt.subplots()
xl = 0
xh = 50
yl = 0
yh = 5
plot_density_grid(
    quadtree,
    fig = fig,
    ax = ax,
    vmin = 1e2,
    vmax = 1e9,
    xl = xl,
    xh = xh,
    yl = yl,
    yh = yh,
)

tau_x_r = range(0, 100, length = 250)
tau_x_z = range(0, 100, length = 250)
tau_x_grid = 1e2 * ones((250, 250))
for (i, r) in enumerate(tau_x_r)
    for (j, z) in enumerate(tau_x_z)
        tau_x = compute_xray_tau(
            smoothed_grid,
            0,
            0,
            r,
            z,
            radiation.xray_luminosity,
            black_hole.Rg,
        )
        tau_x_grid[i, j] = tau_x
    end
end
fig, ax = plt.subplots()
cm = ax.pcolormesh(tau_x_r, tau_x_z, tau_x_grid', norm = LogNorm(1e-3, 10))
plot_streamlines(integrators, fig, ax)
ax.set_xlim(tau_x_r[1], tau_x_r[end])
ax.set_ylim(tau_x_z[1], tau_x_z[end])
plt.colorbar(cm, ax = ax)

@time compute_xray_tau(
    smoothed_grid,
    0,
    0,
    20.0,
    4,
    radiation.xray_luminosity,
    black_hole.Rg,
)

gi = GridIterator(smoothed_grid, 0, 0, 20.0, 4.0)
points = [[0.0, 0.0]]
for point in gi
    push!(points, point)
end
points = hcat(points...)
plt.plot(points[1, :], points[2, :], "o-")



@time compute_uv_tau(
    smoothed_grid,
    0,
    0,
    10.843373493975903,
    8.835341365461847,
    black_hole.Rg,
)

@time compute_uv_tau(
    smoothed_grid,
    1088.7964205397877,
    0,
    6.130665454050429,
    7.999944993615136e-6,
    black_hole.Rg,
)

gi = GridIterator(smoothed_grid, 50, 0, 50, 25)
Qwind.next_intersection!(gi)

r_range = range(0, 100, length = 5)
z_range = range(0, 100, length = 5)
fig, ax = plt.subplots()
for r in r_range
    for z in z_range
        gi = GridIterator(smoothed_grid, 50, 0, r, z)
        println("r $r z $z")
        points = []
        for point in gi
            push!(points, point)
        end
        if length(points) == 0
            continue
        end
        points = hcat(points...)
        println(maximum(points[2, :]))
        ax.plot(points[1, :], points[2, :], "o-")
    end
end

ri = 0
zi = 0
tau_uv_r = range(0, 100, length = 250)
tau_uv_z = range(0, 100, length = 250)
tau_uv_grid = 1e2 * ones((250, 250))
for (i, r) in enumerate(tau_uv_r)
    for (j, z) in enumerate(tau_uv_z)
        println("r $r z $z")
        tau_uv = compute_uv_tau(smoothed_grid, ri, zi, r, z, black_hole.Rg)
        tau_uv_grid[i, j] = tau_uv
    end
end

fig, ax = plt.subplots()
cm = ax.pcolormesh(tau_uv_r, tau_uv_z, tau_uv_grid', norm = LogNorm())
plot_streamlines(integrators, fig, ax)
ax.set_xlim(tau_uv_r[1], tau_uv_r[end])
ax.set_ylim(tau_uv_z[1], tau_uv_z[end])
plt.colorbar(cm, ax = ax)




#plot_grid(quadtree, 10, fig=fig, ax=ax, yh=yh, yl=yl, xl=xl, xh=xh)

#r_range = 10 .^ range(1, log10(3e3), length = 1000);
#z_range = 10 .^ range(-6, 3, length = 1000);
#new_points = hcat([[r, z] for r in r_range for z in z_range]...);
#
#
#
#density_grid_nn =
#    reshape(get_density_points(kdtree, new_points), (length(r_range), length(z_range)))';
#
#fig, ax = plt.subplots()
#cm = ax.pcolormesh(r_range, z_range, density_grid_nn', norm = LogNorm())
#plot_streamlines(integrators, fig, ax)
#plt.colorbar(cm, ax = ax)
#ax.axvline(r_range[end - 20])
#ax.set_xlim(0, 3000)
#ax.set_ylim(0, 1000)
#
#fig, ax = plt.subplots()
#ax.semilogy(z_range, density_grid_nn[end - 20, :])
#
#fig, ax = plt.subplots()
#imgg = imfilter(density_grid_nn, Kernel.gaussian(10))
#cm = ax.pcolormesh(r_range, z_range, imgg', norm = LogNorm())
#plt.colorbar(cm, ax = ax)
#ax.axvline(r_range[end - 20])
##plot_streamlines(integrators, fig, ax)
#ax.set_xlim(0, 3000)
#ax.set_ylim(0, 1000)
#
#fig, ax = plt.subplots()
#ax.semilogy(z_range, imgg[end - 20, :])
#
## interpolte from kdtree
##
#r_range_base = 10 .^ range(1, log10(3e3), length = 250);
#z_range_base = 10 .^ range(-6, 3, length = 250);
#points = hcat([[r, z] for r in r_range_base for z in z_range_base]...);
#densities = get_density_points(kdtree, points);
#
#r_range_final_grid = 10 .^ range(1, log10(3e3), length = 1000);
#z_range_final_grid = 10 .^ range(-6, 3, length = 1000);
#grid_r = r_range' .* ones(length(z_range));
#grid_z = z_range' .* ones(length(r_range));
#
#
#grid_interp = interp.griddata(
#    points',
#    densities,
#    (grid_r, grid_z),
#    method = "cubic",
#    fill_value = 1e2,
#);
#
##
##points = hcat(r, z);
##
#
## grid data
##grid_interp = interp.griddata(points, density, (grid_r, grid_z), method="linear", fill_value=1e2)
#
## interp2d
##int2d = interp.interp2d(r, z, density, fill_value = 1e2, kind="linear")
##grid_interp = int2d(r_range, z_range)
#
## LinearNDInterpolator
##nd_interp = interp.LinearNDInterpolator(points, density, fill_value=1e2)
##grid_interp = nd_interp(grid_r, grid_z)
##
## Rbf
##rbf_interp = interp.Rbf(r, z, density)
##grid_interp = rbf_interp(r_range, z_range)
#
#
##fig, ax = plt.subplots()
##cm = ax.pcolormesh(r_range, z_range, grid_interp', norm = LogNorm())
###plot_streamlines(integrators, fig, ax)
##plt.colorbar(cm, ax=ax)
###ax.scatter(r, z, s=1, color = "white")
##ax.set_xlim(0,2000)
##ax.set_ylim(0,2000)
##grid_interp = rbf_interp(r_range, z_range)
#
#
##fig, ax = plt.subplots()
##cm = ax.pcolormesh(r_range, z_range, grid_interp', norm = LogNorm())
###plot_streamlines(integrators, fig, ax)
##plt.colorbar(cm, ax=ax)
###ax.scatter(r, z, s=1, color = "white")
##ax.set_xlim(0,2000)
##ax.set_ylim(0,2000)
#
#
#


using PyPlot
LogNorm = matplotlib.colors.LogNorm

integrators = iterations_dict[9]["integrators"];
rt = iterations_dict[9]["radiative_transfer"];
kdtree = rt.density_interpolator.kdtree

sg = SmoothedGrid(integrators, kernel_size = 0);

lg = LogGrid(integrators, n_timesteps = 1000)

fig, ax = plt.subplots()
cm = ax.pcolormesh(
    sg.grid.r_range,
    sg.grid.z_range,
    sg.grid.grid',
    norm = LogNorm(1e8, 1e12),
)
plt.colorbar(cm, ax = ax)
plot_streamlines(integrators, fig, ax, color = "black", alpha = 0.5, linestyle = "-")
ax.set_xlim(6, 20)
ax.set_ylim(0, 0.2)
plt.show()

r_range = range(13.5, 14.5, length = 250);
z_range = range(0, 0.02, length = 251);
densities = [get_density(kdtree, r, z) for z in z_range for r in r_range];
densities = reshape(densities, (length(r_range), length(z_range)));
fig, ax = plt.subplots();
cm = ax.pcolormesh(r_range, z_range, densities', norm = LogNorm());
plot_streamlines(integrators, fig, ax, color = "black", alpha = 0.5, linestyle = "o-")
ax.set_xlim(r_range[1], r_range[end])
ax.set_ylim(z_range[1], z_range[end])
plt.colorbar(cm, ax = ax);
plt.show();



# LINE INTERPOLATION

kdtrees = KDTree[]
for integrator in integrators
    r, z, zmax, z0, lw, density =
        get_dense_solution_from_integrator(integrator, 10000)
    kdtree = create_wind_kdtree(r, z, zmax, z0, lw, density, 1e2, 10000)
    push!(kdtrees, kdtree)
end

using ImageFiltering

kdtree = create_wind_kdtree(integrators, 10000)

r_range = 10 .^ range(log10(6), 3, length = 100);
z_range = 10 .^ range(1e-6, log10(10), length = 101);
#densities = [get_density(kdtree, kdtrees, r, z) for z in z_range for r in r_range];
densities = [get_density(kdtree, r, z) for z in z_range for r in r_range];
densities = reshape(densities, (length(r_range), length(z_range)));

#densities_smoothed = imfilter(densities, Kernel.gaussian((7,7), (5,5)));

fig, ax = plt.subplots();
cm = ax.pcolormesh(r_range, z_range, densities', norm = LogNorm());
plot_streamlines(integrators, fig, ax, color = "black", alpha = 0.5, linestyle = "o-")
ax.set_xlim(r_range[1], r_range[end])
ax.set_ylim(z_range[1], z_range[end])
plt.colorbar(cm, ax = ax);
fig.savefig("normal.png", dpi=300)
plt.show();


