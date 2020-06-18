using Qwind
using Test
using Plots
#pyplot()

# create black hole
M = 1e8 * M_SUN
mdot = 0.5
spin = 0.0
black_hole = BlackHole(M, mdot, spin)

# create radiation
f_uv = 0.85
f_x = 0.15
shielding_density = 2e8
r_in = 200.0
radiation = SimpleRadiation(black_hole, f_uv, f_x, shielding_density, r_in)

# create simulation grid
r_min = 0.0
r_max = 10000
z_min = 0.0
z_max = 10000
vacuum_density = 1e2
grid = Grid(r_min, r_max, z_min, z_max)

# initialize streamlines
r_fi = 1600
n_lines = 60
velocity_generator(r_0) = 1e7
density_generator(r_0) = 2e8
streamlines = Streamlines(
    r_in,
    r_fi,
    n_lines,
    velocity_generator,
    density_generator,
    black_hole.M;
    z_0 = 1.0,
    log_spaced = false,
)

#run simulation

for line in streamlines.lines
    parameters = Parameters(line, grid, radiation)
    solver = initialize_solver(line, parameters)
    run_solver!(solver)
end

plt = plot(legend=false)
for line in streamlines.lines
    plot!(line.r, line.z)
end
display(plt)


#tau_x_grid = zeros((2000, 500))
#tau_uv_grid = zeros((2000, 500))
#radiation_r = zeros((2000, 500))
#radiation_z = zeros((2000, 500))
#gravity_r = zeros((2000, 500))
#gravity_z = zeros((2000, 500))
#xi_grid = zeros((2000, 500))
#r_range = range(1, 2000, length=2000) * black_hole.Rg
#z_range = range(1, 10, length=500) * black_hole.Rg
#for (i, r) in enumerate(r_range)
#    for (j, z) in enumerate(z_range)
#        radiation_field_value = compute_disc_radiation_field(radiation, r, z)
#        gravity_value = compute_gravitational_acceleration(r, z, black_hole)
#        gravity_r[i,j] = gravity_value[1]
#        gravity_z[i,j] = gravity_value[2]
#        radiation_r[i,j] = radiation_field_value[1]
#        radiation_z[i,j] = radiation_field_value[2]
#        tau_x_grid[i,j] = compute_xray_tau(radiation, r, z, 2e8, r_in)
#        tau_uv_grid[i,j] = compute_uv_tau(radiation, r, z, 2e8, r_in)
#        xi_grid[i,j] = max(compute_ionization_parameter(radiation, r, z, 2e8, tau_x_grid[i,j]), 1e-20)
#    end
#end
#fig, ax = plt.subplots()
#LogNorm = matplotlib.colors.LogNorm
#cm = plt.pcolormesh(r_range / black_hole.Rg, z_range / black_hole.Rg, transpose(abs.(radiation_z ./ gravity_z)), norm=LogNorm())
#plt.colorbar(cm)
#gcf()

#streamlines = Streamlines(
#    r_in,
#    r_fi,
#    n_lines,
#    velocity_generator,
#    density_generator,
#    black_hole.M,
#    z_0 = black_hole.Rg,
#    log_spaced = false,
#)
#
#line = streamlines.lines[9]
#parameters = Parameters(line, grid, radiation)
#solver = initialize_solver(line, parameters)
