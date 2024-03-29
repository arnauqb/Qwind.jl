module QwindPlotting
export plot_streamlines, plot_density_grid, plot_xray_grid, plot_wind_hull, plot_density_contour, plot_xray_contour, plot_uv_grid
#using PyPlot
#LogNorm = matplotlib.colors.LogNorm
#Normalize = matplotlib.colors.Normalize
#PdfPages = matplotlib.backends.backend_pdf.PdfPages
using PyCall, PyPlot, Qwind
import ConcaveHull
#mpl = pyimport("matplotlib")
#plt = pyimport("matplotlib.pyplot")
pdfpages = pyimport("matplotlib.backends.backend_pdf")
LogNorm = matplotlib.colors.LogNorm
Normalize = matplotlib.colors.Normalize

using ColorSchemes
using ColorSchemes: colorschemes

function plot_streamlines(
    integrators;
    fig = nothing,
    ax = nothing,
    alpha = 1.0,
    xl = nothing,
    xh = nothing,
    yl = nothing,
    yh = nothing,
    linewidth = 1,
    linestyle = "o-",
    colorscheme = nothing,
    color = nothing,
    markersize=nothing
)
    if colorscheme === nothing
        colorscheme_plot = colorschemes[:RdBu_9]
    else
        colorscheme_plot = colorscheme #colorschemes[Symbol(colorscheme)]
    end
    if ax === nothing
        fig, ax = plt.subplots()
    end
    for i = 1:length(integrators)
        if !isassigned(integrators, i)
            continue
        end
        integrator = integrators[i]
        if color === nothing
            color_toplot = get(colorscheme_plot, i / length(integrators))
            color_toplot = [color_toplot.r, color_toplot.g, color_toplot.b]
        else
            color_toplot = color
        end
        ax.plot(
            integrator.p.data[:r],
            integrator.p.data[:z],
            linestyle,
            color = color_toplot,
            alpha = alpha,
            linewidth = linewidth,
            markersize=markersize,
        )
    end
    ax.set_xlabel(L"R [$R_g$]")
    ax.set_ylabel(L"z [$R_g$]")
    if xl !== nothing && xh !== nothing
        ax.set_xlim(xl, xh)
    end
    if yl !== nothing && yh !== nothing
        ax.set_ylim(yl, yh)
    end
    return fig, ax
end

function plot_streamlines(
    integrators::Vector{<:Qwind.Trajectory};
    fig = nothing,
    ax = nothing,
    alpha = 1.0,
    xl = nothing,
    xh = nothing,
    yl = nothing,
    yh = nothing,
    linewidth = 1,
    linestyle = "o-",
    colorscheme = nothing,
    color = nothing,
    markersize=nothing,
)
    if colorscheme === nothing
        colorscheme_plot = colorschemes[:RdBu_9]
    else
        colorscheme_plot = colorscheme #colorschemes[Symbol(colorscheme)]
    end
    if ax === nothing
        fig, ax = plt.subplots()
    end
    for i = 1:length(integrators)
        if !isassigned(integrators, i)
            continue
        end
        integrator = integrators[i]
        if color === nothing
            color_toplot = get(colorscheme_plot, i / length(integrators))
            color_toplot = [color_toplot.r, color_toplot.g, color_toplot.b]
        else
            color_toplot = color
        end
        ax.plot(
            integrator.r,
            integrator.z,
            linestyle,
            color = color_toplot,
            alpha = alpha,
            linewidth = linewidth,
            markersize=markersize,
        )
    end
    ax.set_xlabel(L"R [$R_g$]")
    ax.set_ylabel(L"z [$R_g$]")
    if xl !== nothing && xh !== nothing
        ax.set_xlim(xl, xh)
    end
    if yl !== nothing && yh !== nothing
        ax.set_ylim(yl, yh)
    end
    return fig, ax
end


function plot_density_grid(
    grid::DensityGrid;
    fig = nothing,
    ax = nothing,
    cmap = "viridis",
    vmin = nothing,
    vmax = nothing,
    rmin = 0,
    rmax = 1000,
    zmin = 1e-6,
    zmax = 1000,
    nr = 100,
    nz = 101,
    colorbar=true
)
    if ax === nothing
        fig, ax = plt.subplots()
    end
    r_range = range(rmin, rmax, length=nr)
    z_range = 10 .^ range(log10(zmin), log10(zmax), length=nz)
    r_range_grid = r_range .* ones(nz)'
    z_range_grid = z_range' .* ones(nr)
    density_grid = get_density.(Ref(grid), r_range_grid, z_range_grid)
    cm = ax.pcolormesh(
        r_range,
        z_range,
        density_grid',
        norm = LogNorm(vmin = vmin, vmax = vmax),
        cmap = cmap,
        linewidth=0,
        rasterized=true
    )
    if colorbar
        cb = plt.colorbar(cm, ax = ax)
        cb.set_label(L"Density [cm$^{-3}$]", rotation=-90, labelpad=15)
    end
    return fig, ax, cm
end

function plot_density_contour(
    grid::DensityGrid;
    fig=nothing,
    ax=nothing,
    cmap = "viridis",
    vmin = nothing,
    vmax = nothing,
    xlim = nothing,
    ylim = nothing,
    levels=nothing,
)
    if ax === nothing
        fig, ax = plt.subplots()
    end
    cm = ax.contourf(
        grid.r_range,
        grid.z_range,
        grid.grid',
        norm = LogNorm(),
        cmap = cmap,
        levels=levels,
        shading="auto"
    )
    plt.colorbar(cm, ax = ax)
    if xlim !== nothing
        ax.set_xlim(xlim[1], xlim[2])
    end
    if ylim !== nothing
        ax.set_ylim(ylim[1], ylim[2])
    end
    return fig, ax
end


function plot_xray_grid(
    grid::DensityGrid,
    xray_luminosity,
    Rg;
    fig = nothing,
    ax = nothing,
    cmap = "viridis",
    vmin = nothing,
    vmax = nothing,
    xlim = nothing,
    ylim = nothing,
    rmin = 0,
    rmax = 1000,
    zmin = 1e-6,
    zmax = 1000,
    nr = 100,
    nz = 101,
    zx = 0.0,
    colorbar=true
)
    r_range = range(rmin, rmax, length=nr)
    z_range = 10 .^ range(log10(zmin), log10(zmax), length=nz)
    r_range_grid = r_range .* ones(nz)'
    z_range_grid = z_range' .* ones(nr)
    f(r, z) = compute_tau_xray(grid, Qwind.Boost(), ri=0.0, zi=zx, rf=r, zf=z, xray_luminosity=xray_luminosity, Rg=Rg, mu_nucleon=0.61, mu_electron=1.17)
    ret = f.(r_range_grid, z_range_grid)
    if ax === nothing
        fig, ax = plt.subplots()
    end
    cm = ax.pcolormesh(
        r_range,
        z_range,
        ret',
        norm = LogNorm(vmin = vmin, vmax = vmax),
        cmap = cmap,
        linewidth=0,
        rasterized=true
    )
    if colorbar
        cb = plt.colorbar(cm, ax = ax)
        cb.set_label("X-Ray optical depth", rotation=-90, labelpad=15)
    end
    if xlim !== nothing
        ax.set_xlim(xlim[1], xlim[2])
    end
    if ylim !== nothing
        ax.set_ylim(ylim[1], ylim[2])
    end
    return fig, ax, cm
end


function plot_xray_contour(
    grid::DensityGrid,
    xray_luminosity,
    Rg;
    cmap = "viridis",
    levels = nothing,
    xlim = nothing,
    ylim = nothing,
    rmin = 0,
    rmax = 1000,
    zmin = 1e-6,
    zmax = 1000,
    nr = 100,
    nz = 101,
    zx = 6.0,
    fig=nothing,
    ax=nothing
)
    r_range = range(rmin, rmax, length=nr)
    z_range = 10 .^ range(log10(zmin), log10(zmax), length=nz)
    r_range_grid = r_range .* ones(nz)'
    z_range_grid = z_range' .* ones(nr)
    ret = compute_xray_tau.(Ref(grid), 0.0, zx, r_range_grid, z_range_grid, xray_luminosity, Rg)
    if ax === nothing
        fig, ax = plt.subplots()
    end
    cm = ax.contourf(
        r_range,
        z_range,
        ret',
        norm = LogNorm(),
        levels=levels,
        cmap = cmap,
    )
    plt.colorbar(cm, ax = ax)
    if xlim !== nothing
        ax.set_xlim(xlim[1], xlim[2])
    end
    if ylim !== nothing
        ax.set_ylim(ylim[1], ylim[2])
    end
    return fig, ax
end

function plot_wind_hull(
    hull::ConcaveHull.Hull;
    rmin = 1,
    rmax = 1000,
    zmin = 1e-6,
    zmax = 1000,
    nr = 100,
    nz = 101,
    cmap = "viridis",
    xlim = nothing,
    ylim = nothing,
    fig = nothing,
    ax= nothing,
    alpha=1.0
)
    if ax === nothing
        fig, ax = plt.subplots()
    end
    r_range = 10 .^ range(log10(rmin), log10(rmax), length=nr)
    z_range = 10 .^ range(log10(zmin), log10(zmax), length=nz)
    r_range_grid = r_range .* ones(nz)'
    z_range_grid = z_range' .* ones(nr)

    ret = Qwind.is_point_in_wind.(Ref(hull), r_range_grid, z_range_grid)
    cm = ax.pcolormesh(
        r_range,
        z_range,
        ret',
        cmap = cmap,
        shading="auto",
        linewidth=0,
        alpha=alpha,
        rasterized=true,
        antialiased=true
    )
    if xlim !== nothing
        ax.set_xlim(xlim[1], xlim[2])
    end
    if ylim !== nothing
        ax.set_ylim(ylim[1], ylim[2])
    end
    return fig, ax
end

function plot_uv_grid(
    grid::DensityGrid,
    rd,
    Rg;
    fig = nothing,
    ax = nothing,
    cmap = "viridis",
    vmin = nothing,
    vmax = nothing,
    xlim = nothing,
    ylim = nothing,
    rmin = 0,
    rmax = 1000,
    zmin = 1e-6,
    zmax = 1000,
    nr = 100,
    nz = 101,
    zx = 0.0,
    colorbar=true
)
    r_range = range(rmin, rmax, length=nr)
    z_range = 10 .^ range(log10(zmin), log10(zmax), length=nz)
    r_range_grid = r_range .* ones(nz)'
    z_range_grid = z_range' .* ones(nr)
    ret = compute_uv_tau.(Ref(grid), rd, 0.0, r_range_grid, z_range_grid, Rg)
    if ax === nothing
        fig, ax = plt.subplots()
    end
    cm = ax.pcolormesh(
        r_range,
        z_range,
        ret',
        norm = LogNorm(vmin = vmin, vmax = vmax),
        cmap = cmap,
        linewidth=0,
        rasterized=true
    )
    if colorbar
        cb = plt.colorbar(cm, ax = ax)
        cb.set_label("X-Ray optical depth", rotation=-90, labelpad=15)
    end
    if xlim !== nothing
        ax.set_xlim(xlim[1], xlim[2])
    end
    if ylim !== nothing
        ax.set_ylim(ylim[1], ylim[2])
    end
    return fig, ax, cm
end


end # module

