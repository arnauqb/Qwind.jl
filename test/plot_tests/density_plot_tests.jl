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
    rlim = nothing,
    zlim = nothing,
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
    if rlim === nothing
        r_range_plot = range(minimum(r_range), maximum(r_range), length = nplot)
    else
        r_range_plot = range(rlim[1], rlim[2], length = nplot)
    end
    if zlim === nothing
        z_range_plot = range(minimum(z_range), maximum(z_range), length = nplot)
    else
        z_range_plot = range(zlim[1], zlim[2], length = nplot)
    end
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
    ax[1].scatter(r_range, z_range, color = "white", s = 1)
    ax[1].set_title("Reconstructed")
    ax[2].set_title("Ground truth")
    if rlim === nothing
        ax[1].set_xlim(r_range_plot[1], r_range_plot[end])
    else
        ax[1].set_xlim(rlim)
    end
    if zlim === nothing
        ax[1].set_ylim(z_range_plot[1], z_range_plot[end])
    else
        ax[1].set_ylim(zlim)
    end
    return fig, ax
end

function parabolic_line(
    density_func;
    n = 25,
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
    r_range = range(0, 10, length = n)
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

function random_points(
    density_func;
    n = 25,
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
    r_range = 10 .* rand(n)
    z_range = 10 .* rand(n)
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

function parabolic_lines(
    density_func;
    n = 25,
    log = false,
    vmin = nothing,
    vmax = nothing,
    nplot = 1000,
    atol = 5e-4,
    rtol = 1e-3,
    cell_min_size = 1e-6,
    depth = 5,
    width = 5,
    rlim = nothing,
    zlim = nothing,
)
    r_range_ = range(0, 10, length = n)
    r_range = []
    z_range = []
    parabola(x) = -0.2 * x^2 + 2 * x
    for i = 1:10
        j = 0
        for r in r_range_
            offset = r_range_[end] / length(r_range_) * j
            j += 1
            push!(z_range, parabola(r) + offset)
            push!(r_range, r + i - 1)
        end
    end
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
        rlim = rlim,
        zlim = zlim,
    )
end
