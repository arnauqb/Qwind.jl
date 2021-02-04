using DrWatson
@quickactivate "Qwind"
using RegionTrees, DataFrames, CSV, YAML, Printf
using Qwind
using Profile
using PProf
using PyPlot
LogNorm = matplotlib.colors.LogNorm
using ColorSchemes
using ColorSchemes: colorschemes


function plot_density_grid_nn(
    windkdtree::KDTree;
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
            grid[i, j] = get_density(windkdtree, r, z)
        end
    end
    cm = ax.pcolormesh(r_range, z_range, grid', norm = LogNorm(vmin = vmin, vmax = vmax))
    cbar = plt.colorbar(cm, ax = ax)
    cbar.ax.set_ylabel("Density [cm^-3]", labelpad = 10, rotation = -90)
    ax.set_xlabel("r [Rg]")
    ax.set_ylabel("z [Rg]")
    ax.set_xlim(xl, xh)
    ax.set_ylim(yl, yh)
    return fig, ax
end


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

grid = Grid(config);

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

radiative_transfer = update_radiative_transfer(radiative_transfer, integrators);
qt = radiative_transfer.quadtree;
nn = radiative_transfer.windkdtree;

fig, ax = subplots()
Qwind.plot_streamlines(integrators, fig, ax, colorscheme = "matter");
savefig("results/ryota_meeting/streamlines.png", dpi = 300, bbox_inches = "tight")
show()

fig, ax = subplots()
plot_density_grid_nn(nn, fig = fig, ax = ax, vmin = 1e2, nr = 500, nz = 501);
savefig("results/ryota_meeting/nn_density.png", dpi = 300, bbox_inches = "tight")
show()

fig, ax = subplots()
plot_density_grid_nn(nn, fig = fig, ax = ax, nr = 1000, nz = 1001);
Qwind.plot_streamlines(integrators, fig, ax, color = "white", alpha = 0.5);
savefig("results/ryota_meeting/nn_density_with_lines.png", dpi = 300, bbox_inches = "tight")
show()

fig, ax = subplots()
plot_density_grid(qt, fig = fig, ax = ax, vmin = 1e2, nr = 500, nz = 501);
savefig("results/ryota_meeting/qt_density.png", dpi = 300, bbox_inches = "tight")
show()

fig, ax = subplots()
plot_density_grid(qt, fig = fig, ax = ax, vmin = 1e2, nr = 500, nz = 501);
plot_grid(qt, 7, fig, ax)
savefig("results/ryota_meeting/qt_density_with_grid.png", dpi = 300, bbox_inches = "tight")
show()

fig, ax = subplots()
plot_density_grid(qt, fig = fig, ax = ax, vmin = 1e2, nr = 500, nz = 501);
plot_grid(qt, 7, fig, ax)
Qwind.plot_streamlines(integrators, fig, ax, color = "white", alpha = 0.5);
ax.set_xlim(0, 5000)
ax.set_ylim(0, 5000)
savefig(
    "results/ryota_meeting/qt_density_with_grid_and_lines.png",
    dpi = 300,
    bbox_inches = "tight",
)
show()

fig, ax = subplots()
plot_density_grid(qt, fig = fig, ax = ax, vmin = 1e2, nr = 500, nz = 501);
plot_grid(qt, 10, fig, ax)
tau, coords =
    compute_xray_tau(qt, 1000.0, 250.0, 0.0, 0.0, radiation.xray_luminosity, black_hole.Rg)
coords = hcat(coords...)
ax.plot(coords[1, :], coords[2, :], "o-", color = "blue", linewidth = 1, markersize = 3)
ax.set_xlim(0, 1200)
ax.set_ylim(0, 300)
savefig("results/ryota_meeting/xray_tau_example.png", dpi = 300, bbox_inches = "tight")
show()

