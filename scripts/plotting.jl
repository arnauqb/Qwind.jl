export plot_streamlines, plot_density_grid
#using PyPlot
#LogNorm = matplotlib.colors.LogNorm
#Normalize = matplotlib.colors.Normalize
#PdfPages = matplotlib.backends.backend_pdf.PdfPages
using PyCall, PyPlot
#mpl = pyimport("matplotlib")
#plt = pyimport("matplotlib.pyplot")
pdfpages = pyimport("matplotlib.backends.backend_pdf")
LogNorm = matplotlib.colors.LogNorm
Normalize = matplotlib.colors.Normalize

using ColorSchemes
using ColorSchemes: colorschemes

function plot_streamlines(
    integrators,
    fig = nothing,
    ax = nothing;
    alpha = 1.0,
    xl = nothing,
    xh = nothing,
    yl = nothing,
    yh = nothing,
    linewidth = 1,
    linestyle = "o-",
    colorscheme = nothing,
    color = nothing,
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
            markersize = 1,
            linewidth = linewidth,
        )
    end
    ax.set_xlabel("r [Rg]")
    ax.set_ylabel("z [Rg]")
    if xl !== nothing && xh !== nothing
        ax.set_xlim(xl, xh)
    end
    if yl !== nothing && yh !== nothing
        ax.set_ylim(yl, yh)
    end
    return fig, ax
end


function plot_density_grid(
    grid::InterpolationGrid;
    cmap = "viridis",
    vmin = nothing,
    vmax = nothing,
    xlim = nothing,
    ylim = nothing,
)
    fig, ax = plt.subplots()
    cm = ax.pcolormesh(
        grid.r_range,
        grid.z_range,
        grid.grid',
        norm = LogNorm(vmin = vmin, vmax = vmax),
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
