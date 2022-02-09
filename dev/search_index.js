var documenterSearchIndex = {"docs":
[{"location":"Usage/#Usage","page":"Usage","title":"Usage","text":"","category":"section"},{"location":"Usage/#Configuration","page":"Usage","title":"Configuration","text":"","category":"section"},{"location":"Usage/","page":"Usage","title":"Usage","text":"All the parameters in the code are specified through a parameter file. An example file can be found in configs/config_example.yaml.  All the parameter fields should be self explanatory when familiar with the paper release (TODO: link to arxiv)","category":"page"},{"location":"Usage/","page":"Usage","title":"Usage","text":" black_hole:\n  M: 1e8 # solar masses\n  mdot: 0.5 # Eddington units\n  spin: 0.0 # spin parameter (-1,1)\n\nradiation:\n  relativistic: true \n  f_uv: auto\n  f_x: 0.15\n  n_r: 1000 # number of bins in the disk\n  disk_r_in: 6.0\n  z_xray: 0.0\n  disk_height: 0.0\n  mu_nucleon: 0.61\n  mu_electron: 1.17\n  xray_opacity: boost \n  tau_uv_calculation: no_tau_uv \n  disk_integral_rtol: 1e-3\n  wind_interpolator:\n    update_grid_method: \"average\"\n    vacuum_density: 1e2\n    nz: 500\n    nr: auto \n\ngrid:\n  r_min: 6.0\n  r_max: 50000.0\n  z_min: 0.0\n  z_max: 50000.0\n\ninitial_conditions:\n  mode: CAKIC \n  r_in: 50.0\n  r_fi: 1500.0\n  n_lines: auto\n  log_spaced: true\n  z0: 0.0 # Rg\n  K: 0.03\n  alpha: 0.6\n  use_precalculated: true # use cached results\n\nintegrator:\n  n_iterations: 5\n  atol: 1e-8\n  rtol: 1e-3\n  save_path: \"./example\"","category":"page"},{"location":"Usage/#Running","page":"Usage","title":"Running","text":"","category":"section"},{"location":"Usage/","page":"Usage","title":"Usage","text":"To run the specified parameter file, we do","category":"page"},{"location":"Usage/","page":"Usage","title":"Usage","text":"using Qwind\nmodel = Model(\"configs/config_example.yaml\")\niterations_dict = Dict()\nrun!(model, iterations_dict)","category":"page"},{"location":"Usage/","page":"Usage","title":"Usage","text":"The iterations_dict stores information about each iteration of the density grid, so we know how the wind properties and radiation field are changing throughout iterations. ","category":"page"},{"location":"Usage/#Running-in-parallel","page":"Usage","title":"Running in parallel","text":"","category":"section"},{"location":"Usage/","page":"Usage","title":"Usage","text":"The code automatically runs with all available cores to Julia. The number of cores that Julia can use can be specified as a launch flag, eg","category":"page"},{"location":"Usage/","page":"Usage","title":"Usage","text":"julia -p 6","category":"page"},{"location":"Usage/","page":"Usage","title":"Usage","text":"will launch Julia with 6 cores. Alternatively, we can make us of the Distributed package,","category":"page"},{"location":"Usage/","page":"Usage","title":"Usage","text":"using Distributed\naddprocs(5)","category":"page"},{"location":"Usage/#Reading-the-results","page":"Usage","title":"Reading the results","text":"","category":"section"},{"location":"Usage/","page":"Usage","title":"Usage","text":"In the example above, the results will be stored in ./tests/results.hdf5. A series of functions are built into QWIND to facilitate reading different data of interest.","category":"page"},{"location":"Usage/#Streamlines","page":"Usage","title":"Streamlines","text":"","category":"section"},{"location":"Usage/","page":"Usage","title":"Usage","text":"streamlines = Streamlines(\"example/results.hdf5\");","category":"page"},{"location":"Usage/","page":"Usage","title":"Usage","text":"We can then easily plot them through","category":"page"},{"location":"Usage/","page":"Usage","title":"Usage","text":"using PyPlot\n\nfig, ax = plt.subplots()\nfor sl in sls\n    ax.plot(sl.r, sl.z, linewidth=1)\nend\nax.set_xlim(0,2500)\nax.set_ylim(0,2500)\nax.set_xlabel(L\"$R$ [ $R_g$ ]\")\nax.set_ylabel(L\"$z$ [ $R_g$ ]\")","category":"page"},{"location":"Usage/","page":"Usage","title":"Usage","text":"(Image: Streamlines)","category":"page"},{"location":"Usage/","page":"Usage","title":"Usage","text":"By default, the streamlines of the latest iteration will be loaded. A specific iteration can be checked with","category":"page"},{"location":"Usage/","page":"Usage","title":"Usage","text":"streamlines = Streamlines(\"example/results.hdf5\", 2)","category":"page"},{"location":"Usage/#Density-Grid","page":"Usage","title":"Density Grid","text":"","category":"section"},{"location":"Usage/","page":"Usage","title":"Usage","text":"Similarly,","category":"page"},{"location":"Usage/","page":"Usage","title":"Usage","text":"dgrid = DensityGrid(\"./example/results.hdf5\");\n\nLogNorm = matplotlib.colors.LogNorm\nfig, ax = plt.subplots()\ncm = ax.pcolormesh(dgrid.r_range, dgrid.z_range, dgrid.grid', norm=LogNorm(1e4, 1e8))\nplt.colorbar(cm, ax=ax)\nax.set_xlim(0,2500)\nax.set_ylim(0,2500)\nax.set_xlabel(L\"$R$ [ $R_g$ ]\")\nax.set_ylabel(L\"$z$ [ $R_g$ ]\")","category":"page"},{"location":"Usage/","page":"Usage","title":"Usage","text":"(Image: )","category":"page"},{"location":"Usage/#Velocity-Grid","page":"Usage","title":"Velocity Grid","text":"","category":"section"},{"location":"Usage/","page":"Usage","title":"Usage","text":"Normalize = matplotlib.colors.Normalize\nfig, ax = plt.subplots(1, 3, figsize=(6,2), sharex=true, sharey=true)\nax[1].pcolormesh(vgrid.r_range, vgrid.z_range, vgrid.vr_grid', norm=Normalize(0, 0.5))\nax[2].pcolormesh(vgrid.r_range, vgrid.z_range, vgrid.vphi_grid', norm=Normalize(0, 0.5))\ncm = ax[3].pcolormesh(vgrid.r_range, vgrid.z_range, vgrid.vz_grid', norm=Normalize(0, 0.5))\nplt.colorbar(cm, ax=ax[3])\nax[1].set_ylim(0,2500)\nax[1].set_xlim(0,2500)\nfor i in 1:3\n    ax[i].set_xlabel(L\"$R$ [ $R_g$ ]\")\nend\nax[1].set_ylabel(L\"$z$ [ $R_g$ ]\")\nax[1].set_title(\"Vr\")\nax[2].set_title(\"Vphi\")\nax[3].set_title(\"Vz\")\nplt.subplots_adjust(wspace=0.05, hspace=0.05)","category":"page"},{"location":"Usage/","page":"Usage","title":"Usage","text":"(Image: )","category":"page"},{"location":"example1/#Example-1:-Mdot-dependence-with-radius","page":"Example 1: Mdot dependence with radius","title":"Example 1: Mdot dependence with radius","text":"","category":"section"},{"location":"example1/","page":"Example 1: Mdot dependence with radius","title":"Example 1: Mdot dependence with radius","text":"This example illustrates how one can change the normalised accretion rate as a function of radius.","category":"page"},{"location":"example1/","page":"Example 1: Mdot dependence with radius","title":"Example 1: Mdot dependence with radius","text":"Let us load our example model","category":"page"},{"location":"example1/","page":"Example 1: Mdot dependence with radius","title":"Example 1: Mdot dependence with radius","text":"using Qwind\nmodel = Model(\"configs/config_example.yaml\")","category":"page"},{"location":"example1/","page":"Example 1: Mdot dependence with radius","title":"Example 1: Mdot dependence with radius","text":"since we specified mdot=0.5, the mdot grid is flat at 0.5:","category":"page"},{"location":"example1/","page":"Example 1: Mdot dependence with radius","title":"Example 1: Mdot dependence with radius","text":"using PyPlot\nplt.plot(model.rad.disk_grid, model.rad.mdot_grid)","category":"page"},{"location":"example1/","page":"Example 1: Mdot dependence with radius","title":"Example 1: Mdot dependence with radius","text":"Let us now run the model for a couple of iterations","category":"page"},{"location":"example1/","page":"Example 1: Mdot dependence with radius","title":"Example 1: Mdot dependence with radius","text":"iterations_dict = Dict()\nrun!(model, iterations_dict, start_it=1, n_iterations=2)","category":"page"},{"location":"example1/","page":"Example 1: Mdot dependence with radius","title":"Example 1: Mdot dependence with radius","text":"We can plot the last iteration,","category":"page"},{"location":"example1/","page":"Example 1: Mdot dependence with radius","title":"Example 1: Mdot dependence with radius","text":"streamlines = iterations_dict[2][\"streamlines\"]\nfig, ax = plt.subplots()\nfor sl in streamlines\n    ax.plot(sl.r, sl.z)\nend","category":"page"},{"location":"example1/","page":"Example 1: Mdot dependence with radius","title":"Example 1: Mdot dependence with radius","text":"At this point, imagine we would like to change our mdot grid to have a different value with radius. Let's suppose we want it to go from 0.2 to 0.5, we can do:","category":"page"},{"location":"example1/","page":"Example 1: Mdot dependence with radius","title":"Example 1: Mdot dependence with radius","text":"n_disk = length(model.rad.disk_grid)\nnew_mdot_grid = range(0.2, 0.5, length=n_disk)\nmodel.rad.mdot_grid .= new_mdot_grid","category":"page"},{"location":"example1/","page":"Example 1: Mdot dependence with radius","title":"Example 1: Mdot dependence with radius","text":"And check it has updated properly","category":"page"},{"location":"example1/","page":"Example 1: Mdot dependence with radius","title":"Example 1: Mdot dependence with radius","text":"plt.plot(model.rad.disk_grid, model.rad.mdot_grid)","category":"page"},{"location":"example1/","page":"Example 1: Mdot dependence with radius","title":"Example 1: Mdot dependence with radius","text":"Now we can run one more iteration","category":"page"},{"location":"example1/","page":"Example 1: Mdot dependence with radius","title":"Example 1: Mdot dependence with radius","text":"run!(model, iterations_dict, start_it=3, n_iterations=1)","category":"page"},{"location":"example1/","page":"Example 1: Mdot dependence with radius","title":"Example 1: Mdot dependence with radius","text":"and check the new results","category":"page"},{"location":"example1/","page":"Example 1: Mdot dependence with radius","title":"Example 1: Mdot dependence with radius","text":"streamlines = iterations_dict[3][\"streamlines\"]\nfig, ax = plt.subplots()\nfor sl in streamlines\n    ax.plot(sl.r, sl.z)\nend","category":"page"},{"location":"Setup/#Setup","page":"Setup","title":"Setup","text":"","category":"section"},{"location":"Setup/","page":"Setup","title":"Setup","text":"The code is written in the Julia programming language. If you are not familiar with it, don't worry! To run the code at a higher level you will need very little Julia knowledge, although I would encourage you to use it for your scientific computing projects. Refer to the official web-page for installation instructions.","category":"page"},{"location":"Setup/","page":"Setup","title":"Setup","text":"The easiest way to install QWIND is through the Julia package index.","category":"page"},{"location":"Setup/","page":"Setup","title":"Setup","text":"using Pkg\nPkg.add(\"Qwind\")","category":"page"},{"location":"Setup/","page":"Setup","title":"Setup","text":"This will automatically install any dependencies you may need. ","category":"page"},{"location":"change_fm/#Using-a-different-force-multiplier","page":"Using a different force multiplier","title":"Using a different force multiplier","text":"","category":"section"},{"location":"change_fm/","page":"Using a different force multiplier","title":"Using a different force multiplier","text":"In this section we explain how to change the algorithm that computes the force multiplier. The computation needs to take 2 arguments: effective optical depth t and ionisation parameter xi.","category":"page"},{"location":"change_fm/","page":"Using a different force multiplier","title":"Using a different force multiplier","text":"We first need to import Qwind and the compute_force_multiplier function","category":"page"},{"location":"change_fm/","page":"Using a different force multiplier","title":"Using a different force multiplier","text":"using Qwind\nimport Qwind: compute_force_mulitplier","category":"page"},{"location":"change_fm/","page":"Using a different force multiplier","title":"Using a different force multiplier","text":"We now need to define a new type to dispatch our calculation of the force multiplier to.","category":"page"},{"location":"change_fm/","page":"Using a different force multiplier","title":"Using a different force multiplier","text":"struct FMNewFlag <: FMInterpolationFlag end","category":"page"},{"location":"change_fm/","page":"Using a different force multiplier","title":"Using a different force multiplier","text":"and then we can define our new method:","category":"page"},{"location":"change_fm/","page":"Using a different force multiplier","title":"Using a different force multiplier","text":"function compute_force_multiplier(t, ionization_parameter, mode::FMNewFlag)\n    return 1.0\nend","category":"page"},{"location":"change_fm/","page":"Using a different force multiplier","title":"Using a different force multiplier","text":"Now we need to specify the new mode in the parameters that are passed to the model. For now the way to do it is:","category":"page"},{"location":"change_fm/","page":"Using a different force multiplier","title":"Using a different force multiplier","text":"parameters = Parameters(\"../configs/config_example.yaml\")\nparameters = change_parameter(parameters, :fm_interp_method_flag, FMNewFlag()) # For now new modes cannot be specified in the .yaml files","category":"page"},{"location":"change_fm/","page":"Using a different force multiplier","title":"Using a different force multiplier","text":"and then we just initialise and run the model","category":"page"},{"location":"change_fm/","page":"Using a different force multiplier","title":"Using a different force multiplier","text":"model = Model(parameters)\niterations_dict = run!(model)","category":"page"},{"location":"change_fm/","page":"Using a different force multiplier","title":"Using a different force multiplier","text":"And to make sure it worked we can plot the force multiplier for one streamline:","category":"page"},{"location":"change_fm/","page":"Using a different force multiplier","title":"Using a different force multiplier","text":"using Plots\nintegrator = iterations_dict[2][\"integrators\"][10]\nfm = integrator.p.data[:fm]\nplot(fm)","category":"page"},{"location":"#QWIND","page":"QWIND","title":"QWIND","text":"","category":"section"},{"location":"","page":"QWIND","title":"QWIND","text":"QWIND is a computer code that models radiation-driven winds originating from accretion discs.","category":"page"},{"location":"","page":"QWIND","title":"QWIND","text":"The wind is modelled as a set of gas trajectories initialising from the surface of the accretion disc. The trajectories are computed by solving the equation of motion given by the radiation and gravitational force. To compute the radiation force, the radiation field is simulated across the wind, taking into account the disc as an extended source and a point-like source at the centre which aims to represent the X-ray corona.","category":"page"}]
}
