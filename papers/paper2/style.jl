using ColorSchemes
using Colors
using PyPlot
LogNorm = matplotlib.colors.LogNorm
plt.style.use("science")


cmap_hex = ["f94144","f3722c","f8961e","f9c74f","90be6d","43aa8b","577590"]
#cmap_hex = ["264653","2a9d8f","e9c46a","f4a261","e76f51"]
cmap_hex =["006466","065a60","0b525b","144552","1b3a4b","212f45","272640","312244","3e1f47","4d194d"]

cmap_parsed = [parse(Colorant, "#$(color)") for color in cmap_hex]
gradient = ColorScheme(cmap_parsed, "paper scheme")
cmap = ColorMap("paper", cmap_parsed)

