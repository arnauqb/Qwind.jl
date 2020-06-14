using Qwind
using Test
using PyPlot

# create black hole
M = 1e8 * M_SUN
mdot = 0.5
spin = 0.0
black_hole = BlackHole(M, mdot, spin)

# create radiation
f_uv = 0.85
f_x = 0.15
uv_fractions = []
radiation = Radiation(black_hole, f_uv, f_x, uv_fractions)

# create simulation grid
r_min = 0.0
r_max = 2000 * Rg(black_hole)
z_min = 0.0
z_max = 2000 * Rg(black_hole)
vacuum_density = 1e2
grid = Grid(r_min, r_max, z_min, z_max, vacuum_density)

# initialize streamlines
r_in = 50 * Rg(black_hole)
lines_r_range = 10 .^ range(log10(r_in), log10(1500.0), length=20.0)

