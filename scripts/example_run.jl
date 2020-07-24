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

# radiative transfer model 
shielding_density = 2e8
rin = 50.0
radiative_transfer = RERadiativeTransfer(radiation, shielding_density, rin, black_hole.Rg)

# create simulation grid
r_min = 0.0
r_max = 10000
z_min = 1.0
z_max = 10000
grid = Grid(r_min, r_max, z_min, z_max)

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
    grid,
    initial_conditions,
    atol = 1e-8,
    rtol = 1e-3,
)

##
# run integrators
run_integrators!(integrators)

# plot lines
#fig, ax = plot_streamlines(integrators, xh = 10000, yh = 2000)
#gcf()

adaptive_mesh = AdaptiveMesh(radiation, integrators)


##
# plot density 



##
#integrators_dict = Dict()

#for i in 1:10
#    global integrators
#    global integrators_dict
#    # construct mesh
#    adaptive_mesh = AdaptiveMesh(radiation, integrators)
#
#    # initialize itnegrators
#    integrators = initialize_integrators(
#        adaptive_mesh,
#        grid,
#        initial_conditions,
#        atol = 1e-8,
#        rtol = 1e-3,
#    )
#    # run integrators
#    run_integrators!(integrators)
#    integrators_dict[i] = integrators
#end
