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
r_fi = 1500 * Rg(black_hole)
n_lines = 30
velocity_generator(r_0) = 5e7
density_generator(r_0) = 2e8
streamlines = Streamlines(
    r_in,
    r_fi,
    n_lines,
    velocity_generator,
    density_generator,
    black_hole.M;
    z_0 = 0.0,
    log_spaced = true,
)


