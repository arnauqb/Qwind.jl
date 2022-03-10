# Running a parameter grid

The main goal of Qwind is to be able to run fast parameter searches. To make it easy, the user can specify ranges to vary in the parameter file. 

Let's illustrate this with an example changing BH mass and mass accretion rate. The parameter file we use (located at `configs/parameter_scan_example.yaml` ) is



```yaml
system:
  qwind_path: /path/to/Qwind.jl
  n_cpus: 16
  partition: cosma
  account: durham
  max_time: 72:00:00
  job_name: qwind

black_hole:
  M: "@vary grid 1e8,1e9,1e10" # take the values specified
  mdot: "@vary log 0.1 0.5 5" # from 0.1 to 0.5, 5 values log spaced
  spin: 0.0

radiation:
  relativistic: true
  f_uv: 0.85
  f_x: 0.15
  n_r: 2000
  disk_r_in: 6.0
  z_xray: 0.0
  disk_height: 0.0
  uv_opacity: thomson
  xray_opacity: boost 
  xray_scattering: false
  tau_uv_calculation: no_tau_uv
  vacuum_density: 1e2
  update_grid_method: "average"
  grid_nz: 250
  grid_nr: auto
  mu_nucleon: 0.61

initial_conditions:
  mode: cak
  r_in: 20.0
  r_fi: 1500.0
  n_trajs: 50
  z0: 0.0 
  K: 0.03
  alpha: 0.6
  use_precalculated: false

integrator:
  n_iterations: 5
  save_path: "./parameter_scan_example"
  integrator_r_min: 6.0
  integrator_r_max: 10000.0
  integrator_z_min: 0.0
  integrator_z_max: 10000.0

numerical_tolerances:
  disk_integral_atol: 0
  disk_integral_rtol: 1e-3
  disk_integral_maxevals: 10000
  integrator_atol: 1e-8
  integrator_rtol: 1e-3
  scattering_atol: 0
  scattering_rtol: 1e-3
```

The first config blog is the `system` configuration block, which specifies the path to our Qwind installation, the number of cpus to use for each model, the partition and account of our HPC cluster, the maximum time of the job, and the job name. Note that this only works for a SLURM scheduler.

We specify the BH mass and mass accretion rate to vary in a certain range. For the BH mass we have specified a grid of values, and for mdot we specify a log range. Another possible option could have been `@vary linear 0.1 0.5 5` which would have varied the value linearly rather than in log-space. 

We can now create a folder hierarchy structure with the variations like this

```julia
using Qwind
create_models_folders("configs/parameter_scan_example.yaml")
```

This will create a folder called `parameter_scan_example` with the following structure

```
├── all_configs.yaml
|
├── model_001
│   ├── config.yaml
│   ├── submit.sh
├── model_002
│   ├── config.yaml
│   ├── submit.sh
...
├── run_model.jl
├── stdout
├── submit_all.sh
└── all_configs.yaml
```

In total we can see 15 models, as we would expect from 3 x 5. Each model has its own `config.yaml` plus a submit script for our system. We can then submit the scripts one by one, or submit them all at once with `bash submit_all.sh`. The logs of the jobs will be stored in the `stdout` folder.

