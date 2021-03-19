using ColorSchemes, Colors, PyPlot
LogNorm = matplotlib.colors.LogNorm
Normalize = matplotlib.colors.Normalize
plt.style.use("./mnras.mplstyle")


cmap_hex = ["f94144","f3722c","f8961e","f9c74f","90be6d","43aa8b","577590"]
#cmap_hex = ["264653","2a9d8f","e9c46a","f4a261","e76f51"]
cmap_hex =["006466","065a60","0b525b","144552","1b3a4b","212f45","272640","312244","3e1f47","4d194d"]

cmap_parsed = [parse(Colorant, "#$(color)") for color in cmap_hex]
gradient = ColorScheme(cmap_parsed, "paper scheme")
cmap = ColorMap("paper", cmap_parsed)

function set_size(width="mnras"; fraction=1, subplots=(1, 1), double=false)
    """Set figure dimensions to avoid scaling in LaTeX.

    Parameters
    ----------
    width: float or string
            Document width in points, or string of predined document type
    fraction: float, optional
            Fraction of the width which you wish the figure to occupy
    subplots: array-like, optional
            The number of rows and columns of subplots.
    Returns
    -------
    fig_dim: tuple
            Dimensions of figure in inches
    """
    if width == "mnras"
        if double
            width_pt = 504
        else
            width_pt = 240
        end
    else
        width_pt = width
    end

    # Width of figure (in pts)
    fig_width_pt = width_pt * fraction
    # Convert from pt to inches
    inches_per_pt = 1 / 72.27

    # Golden ratio to set aesthetic figure height
    # https://disq.us/p/2940ij3
    golden_ratio = (5^0.5 - 1) / 2

    # Figure width in inches
    fig_width_in = fig_width_pt * inches_per_pt
    # Figure height in inches
    fig_height_in = fig_width_in * golden_ratio * (subplots[1] / subplots[2])

    return (fig_width_in, fig_height_in)
end
