using PyPlot
using RegionTrees
LogNorm = matplotlib.colors.LogNorm
Normalize = matplotlib.colors.Normalize
using ColorSchemes
using ColorSchemes: colorschemes

using Qwind

function plot_density_grid(
    r_range,
    z_range,
    density_grid;
    fig = nothing,
    ax = nothing,
    colorbar = false,
    vmin = nothing,
    vmax = nothing,
    log = true,
)
    if ax === nothing
        fig, ax = plt.subplots()
    end
    if log
        norm = LogNorm(vmin = vmin, vmax = vmax)
    else
        norm = Normalize(vmin = vmin, vmax = vmax)
    end
    cm = ax.pcolormesh(r_range, z_range, density_grid', norm = norm)
    if colorbar
        plt.colorbar(cm, ax = ax)
    end
    ax.set_xlabel("r [Rg]")
    ax.set_ylabel("z [Rg]")
    return fig, ax
end

function make_density_grid(r_range, z_range, density_getter;)
    grid_den = zeros(length(r_range), length(z_range))
    for (i, r) in enumerate(r_range[1:(end - 1)])
        for (j, z) in enumerate(z_range[1:(end - 1)])
            grid_den[i, j] = density_getter(r, z)[1]
        end
    end
    return grid_den
end

function make_density_grid(r_range, z_range, quadtree::Cell)
    density_getter = (r, z) -> [get_density(quadtree, r, z)]
    return make_density_grid(r_range, z_range, density_getter)
end

function compare_density_grids(
    density_func,
    r_range,
    z_range;
    log = false,
    vmin = nothing,
    vmax = nothing,
    nplot = 1000,
    atol = 5e-4,
    rtol = 1e-3,
    cell_min_size = 1e-6,
    depth = 5,
    width = 5,
)
    width_range = width .* ones(length(r_range))
    quadtree, bh = create_test_quadtree(
        density_func,
        r_range = r_range,
        z_range = z_range,
        atol = atol,
        rtol = rtol,
        width_range = width_range,
    )
    vmin = vmin
    vmax = vmax
    fig, ax = plt.subplots(1, 2)
    r_range_plot = range(0, 10, length = nplot)
    z_range_plot = range(0, 10, length = nplot)
    grid_quadtree = make_density_grid(r_range_plot, z_range_plot, quadtree)
    func_quadtree = make_density_grid(r_range_plot, z_range_plot, density_func)
    plot_density_grid(
        r_range_plot,
        z_range_plot,
        grid_quadtree,
        fig = fig,
        ax = ax[1],
        vmin = vmin,
        vmax = vmax,
        log = log,
    )
    plot_density_grid(
        r_range_plot,
        z_range_plot,
        func_quadtree,
        fig = fig,
        ax = ax[2],
        colorbar = true,
        vmin = vmin,
        vmax = vmax,
        log = log,
    )
    plot_grid(quadtree, depth, fig, ax[1])
    ax[1].plot(r_range, z_range, "o-", color = "white", linestyle = "-")
    ax[1].set_title("Reconstructed")
    ax[2].set_title("Ground truth")
    ax[1].set_xlim(r_range_plot[1], r_range_plot[end])
    ax[2].set_xlim(z_range_plot[1], z_range_plot[end])
    return fig, ax
end

function parabolic_line(
    density_func;
    log = false,
    vmin = nothing,
    vmax = nothing,
    nplot = 1000,
    atol = 5e-4,
    rtol = 1e-3,
    cell_min_size = 1e-6,
    depth = 5,
    width = 5,
)
    r_range = range(0, 10, length = 10)
    parabola(x) = -0.2 * x^2 + 2 * x
    z_range = parabola.(r_range)
    return compare_density_grids(
        density_func,
        r_range,
        z_range,
        log = log,
        vmin = vmin,
        vmax = vmax,
        nplot = nplot,
        atol = atol,
        rtol = rtol,
        cell_min_size = cell_min_size,
        depth = depth,
        width = width,
    )
end
