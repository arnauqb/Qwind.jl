using PyPlot
LogNorm = matplotlib.colors.LogNorm
using ColorSchemes
using ColorSchemes: colorschemes
export plot_grid,
    plot_cell,
    plot_streamlines,
    plot_density,
    plot_taux_grid,
    plot_density_grid_nn,
    plot_density_grid

function plot_grid(
    quadtree::Cell,
    depth;
    fig = nothing,
    ax = nothing,
    xl = nothing,
    xh = nothing,
    yl = nothing,
    yh = nothing,
)
    if ax === nothing
        fig, ax = plt.subplots()
    end
    children_list = Any[]
    currentchildren = Any[]
    for kid in quadtree.children
        push!(currentchildren, kid)
        push!(children_list, kid)
    end
    currentchildren_aux = Any[]
    dp = depth
    while dp > 0
        for kid in currentchildren
            if kid.children === nothing
                continue
            end
            for kid2 in kid.children
                push!(children_list, kid2)
                push!(currentchildren_aux, kid2)
            end
        end
        currentchildren = copy(currentchildren_aux)
        currentchildren_aux = []
        dp -= 1
    end
    for cell in children_list
        v = hcat(collect(vertices(cell.boundary))...)
        x = v[1, [1, 2, 4, 3, 1]]
        y = v[2, [1, 2, 4, 3, 1]]
        ax.plot(x, y, color = "black", alpha = 0.1)
    end
    ax.set_xlabel("r [Rg]")
    ax.set_ylabel("z [Rg]")
    if xl !== nothing
        ax.set_xlim(xl, xh)
    end
    if yl !== nothing
        ax.set_ylim(yl, yh)
    end
    return fig, ax
end

function plot_streamlines(
    integrators,
    fig = nothing,
    ax = nothing;
    alpha = 1.0,
    xl = 0,
    xh = 2000,
    yl = 0,
    yh = 2000,
    linewidth=1,
    colorscheme = nothing,
    color = nothing, 
)
    if colorscheme === nothing
        colorscheme_plot = colorschemes[:RdBu_9]
    else
        colorscheme_plot = colorschemes[Symbol(colorscheme)]
    end
    if ax === nothing
        fig, ax = plt.subplots()
    end
    for (i, integrator) in enumerate(integrators)
        if color === nothing
            color_toplot = get(colorscheme_plot, i / length(integrators))
            color_toplot = [color_toplot.r, color_toplot.g, color_toplot.b]
        else
            color_toplot = color
        end
        ax.plot(
            integrator.p.data.r,
            integrator.p.data.z,
            "o-",
            color = color_toplot,
            alpha = alpha,
            markersize = 1,
            linewidth=linewidth,
        )
    end
    ax.set_xlabel("r [Rg]")
    ax.set_ylabel("z [Rg]")
    ax.set_xlim(xl, xh)
    ax.set_ylim(yl, yh)
    return fig, ax
end

function plot_density(
    quadtree::Cell;
    xl = 0,
    xh = 5000,
    yl = 0,
    yh = 5000,
    grid = true,
    depth = 4,
    nr = 500,
    nz = 501,
    vmin = nothing,
    vmax = nothing,
    fig = nothing,
    ax = nothing,
)
    if ax === nothing
        fig, ax = plt.subplots()
    end
    grid_den = zeros(nr, nz)
    r_range = range(xl, stop = xh, length = nr)
    z_range = range(yl, stop = yh, length = nz)
    for (i, r) in enumerate(r_range[1:(end - 1)])
        for (j, z) in enumerate(z_range[1:(end - 1)])
            grid_den[i, j] = get_density(quadtree, r, z)
        end
    end
    cm = ax.pcolormesh(r_range, z_range, grid_den', norm = LogNorm())
    plt.colorbar(cm, ax = ax)
    ax.set_xlabel("r [Rg]")
    ax.set_ylabel("z [Rg]")
    ax.set_xlim(xl, xh)
    ax.set_ylim(yl, yh)
    return fig, ax
end

#
#function plot_density(
#    windtree::WindTree;
#    xl = 0,
#    xh = 5000,
#    yl = 0,
#    yh = 5000,
#    grid = true,
#    depth = 4,
#    nr = 500,
#    nz = 501,
#    vmin = nothing,
#    vmax = nothing,
#    fig = nothing,
#    ax = nothing,
#)
#    if ax === nothing
#        fig, ax = plt.subplots()
#    end
#    grid_den = zeros(nr, nz)
#    r_range = range(xl, stop = xh, length = nr)
#    z_range = range(yl, stop = yh, length = nz)
#    for (i, r) in enumerate(r_range[1:end-1])
#        for (j, z) in enumerate(z_range[1:end-1])
#            grid_den[i, j] = get_density_from_tree(windtree, r, z)
#        end
#    end
#    cm = ax.pcolormesh(r_range, z_range, grid_den', norm = LogNorm())
#    plt.colorbar(cm, ax = ax)
#    ax.set_xlabel("r [Rg]")
#    ax.set_ylabel("z [Rg]")
#    ax.set_xlim(xl, xh)
#    ax.set_ylim(yl, yh)
#    return fig, ax
#end
#
#function plot_taux_grid(
#    quadtree,
#    xray_luminosity,
#    Rg;
#    z_0 = 0.0,
#    xl = 0,
#    xh = 5000,
#    yl = 0,
#    yh = 5000,
#    grid = true,
#    depth = 4,
#    nr = 500,
#    nz = 501,
#    vmin = nothing,
#    vmax = nothing,
#    plt = nothing,
#)
#    if plt === nothing
#        plt = plot()
#    end
#    grid = zeros(Float64, nr, nz)
#    r_range = range(xl, stop = xh, length = nr)
#    z_range = range(yl, stop = yh, length = nz)
#    for (i, r) in enumerate(r_range[1:end-1])
#        for (j, z) in enumerate(z_range[1:end-1])
#            grid[i, j] =
#                compute_xray_tau(quadtree::Cell, r, z, z_0, xray_luminosity, Rg)
#        end
#    end
#    heatmap!(plt, r_range, z_range, log10.(grid)')
#    xlabel!("r [Rg]")
#    ylabel!("z [Rg]")
#    xlims!((xl, xh))
#    ylims!(yl, yh)
#    return plt
#end
#
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
    plt.colorbar(cm, ax = ax)
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
