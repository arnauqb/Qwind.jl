# imports
##
using DrWatson
@quickactivate "Qwind"
using Qwind

using PyPlot
LogNorm = matplotlib.colors.LogNorm
using RegionTrees


# create black hole
##
M = 1e8 * M_SUN
mdot = 0.5
spin = 0.0
black_hole = BlackHole(M, mdot, spin)

# create radiation model
##
fuv = 0.85
fx = 0.15
radiation = RERadiation(black_hole, fuv, fx)
#radiation = QsosedRadiation(black_hole, 300, fx)

# radiative transfer model
shielding_density = 2e8
rin = 50.0
radiative_transfer = RERadiativeTransfer(radiation, shielding_density, rin, black_hole.Rg)

# create simulation grid
r_min = 0.0
r_max = 10000
z_min = 1.0
z_max = 10000
sim_grid = Grid(r_min, r_max, z_min, z_max)

# initial conditions
z0 = 1.0
n0 = 2e8
v0 = 5e7 / C
rfi = 1600
nlines = 100
logspaced = true
initial_conditions = UniformIC(rin, rfi, nlines, z0, n0, v0, logspaced)


##
# initialize itnegrators
integrators = initialize_integrators(
    radiative_transfer,
    sim_grid,
    initial_conditions,
    atol = 1e-8,
    rtol = 1e-3,
)

iterations = Dict()
iterations[0] = integrators

iterations_dict = Dict(0 => integrators)
##
# run integrators
run_integrators!(integrators)

# plot lines
#fig, ax = plot_streamlines(integrators, xh = 10000, yh = 2000)
#gcf()

# iteration 2
radiation = QsosedRadiation(black_hole, 300, fx)

tau_max_std = 0.01

adaptive_mesh = AdaptiveMesh(radiation, integrators, tau_max_std=tau_max_std)

initial_conditions = CAKIC(radiation, black_hole, rin, rfi, nlines, z0, 0.03, 0.6, 0.5, logspaced)

integrators = initialize_integrators(adaptive_mesh, sim_grid, initial_conditions, atol=1e-8, rtol=1e-3)

iterations[1] = integrators


run_integrators!(integrators)

adaptive_mesh = AdaptiveMesh(radiation, integrators, tau_max_std=tau_max_std)

integrators = initialize_integrators(adaptive_mesh, sim_grid, initial_conditions, atol=1e-8, rtol=1e-3)

iterations[2] = integrators

run_integrators!(integrators)

adaptive_mesh = AdaptiveMesh(radiation, integrators, tau_max_std=tau_max_std)

# iteration 3

#end
# start iterations
#n_iterations = 3
#
#for i in 1:n_iterations
#    global adaptive_mesh
#    global iterations
#    global tau_max_std
#    integrators = initialize_integrators(adaptive_mesh, sim_grid, initial_conditions, atol=1e-8, rtol=1e-3)
#    iterations[i] = integrators
#    run_integrators!(integrators)
#    adaptive_mesh = AdaptiveMesh(radiation, integrators, tau_max_std=tau_max_std)
#end
