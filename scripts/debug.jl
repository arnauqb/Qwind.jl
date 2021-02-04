using DrWatson
@quickactivate "Qwind"
using Qwind
using Profile
using PProf, CSV, DataFrames

#model = Model("scripts/config.yaml")
#iterations_dict = Dict()
#run!(model, iterations_dict)

lines_range, lines_widths = get_initial_radii_and_linewidths(model.ic, model.bh.Rg);

rt = iterations_dict[4]["radiative_transfer"]

line_id = 15;

integrator = initialize_integrator(
    rt,
    model.wind_grid,
    model.ic,
    lines_range[15],
    lines_widths[15],
    atol = 1e-6,
    rtol = 1e-3,
    line_id = 15,
    save_path = "tests",
);

run_integrator!(integrator)

plot(integrator.p.data[:r], integrator.p.data[:z])

semilogy(integrator.p.data[:vr])

r = integrator.p.data[:r][end]
z = integrator.p.data[:z][end]


compute_gravitational_acceleration(r, z)

compute_disc_radiation_field(model.rt, r, z)

integrator.p.l0^2 / r^3
