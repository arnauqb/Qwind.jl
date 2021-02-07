using DrWatson
@quickactivate "Qwind"
using Qwind
using Profile
using PProf, CSV, DataFrames

model = Model("paper2/plot_nn.yaml")
iterations_dict = Dict()
#run!(model, iterations_dict)

lines_range, lines_widths = get_initial_radii_and_linewidths(model.ic, model.bh.Rg);

line_id = 75;
integrator = initialize_integrator(
    model.rt,
    model.wind_grid,
    model.ic,
    lines_range[line_id],
    lines_widths[line_id],
    atol = 1e-8,
    rtol = 1e-3,
    line_id = line_id,
    save_path = "tests",
);
run_integrator!(integrator);
loglog(integrator.p.data[:z], integrator.p.data[:n])

#loglog(integrator.p.data[:r], integrator.p.data[:z])

semilogy(integrator.p.data[:vr])

r = integrator.p.data[:r][end]
z = integrator.p.data[:z][end]


compute_gravitational_acceleration(r, z)

compute_disc_radiation_field(model.rt, r, z)

integrator.p.l0^2 / r^3
