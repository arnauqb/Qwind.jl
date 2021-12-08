using Distributed
@everywhere using DrWatson
@everywhere @quickactivate "Qwind"
@everywhere using Qwind

using YAML, HDF5, CSV, DataFrames, PyPlot
include("scripts/plotting.jl")
using Profile, PProf, TimerOutputs, BenchmarkTools, ProgressMeter
LogNorm = matplotlib.colors.LogNorm
Normalize = matplotlib.colors.Normalize

function get_model(config)
    model = Model(config)
    try
        mv(model.parameters.save_path, "backup", force = true)
    catch
    end
    model = Model(config)
    iterations_dict = Dict()
    return model, iterations_dict
end
model, iterations_dict = get_model("./configs/proga.yaml");
run!(model, iterations_dict, parallel = true)

rr = range(60, 1500, length=500);
n0s = getn0.(Ref(model), rr);

fig, ax = plt.subplots()
ax.loglog(rr, n0s)


Qwind.compute_proga_density(model.bh, 100)

getn0(model, 100)

integrators = iterations_dict[5]["integrators"];
times, merging = Qwind.get_intersection_times(integrators);



data = CSV.read("/home/arnau/code/qwind_experiments/data/james_rho.csv", DataFrame);
data[!, :r] = data.r / model.bh.Rg;
data[!, :z] = data.z / model.bh.Rg;
data[!, :rho] = data.rho / (model.parameters.mu_nucleon * M_P);
points = hcat(log10.(data.r), log10.(data.z));
interp = Qwind.scipy_interpolate.LinearNDInterpolator(points, log10.(data.rho), fill_value=2);
r_range = 10 .^ range(log10(minimum(data[!, :r])), log10(maximum(data[!, :r])), length=250);
z_range = 10 .^ range(log10(minimum(data[!, :z])), log10(maximum(data[!, :z])), length=251);
values = zeros(length(r_range), length(z_range));
for (i, r) in enumerate(r_range)
    for (j, z) in enumerate(z_range)
        values[i, j] = max(1e2, 10 .^ interp(log10(r), log10(z))[1])
    end
end
density_grid = DensityGrid(r_range, z_range, values);

r_range = 10 .^ range(log10(60), log10(1500), length=250);
density_profile = []
for r in r_range
    value = 1e2
    z = 0
    while value == 1e2
        value = get_density(density_grid, r, z)
        z += 1e-2
        if z > 50
            break
        end
    end
    push!(density_profile, value)
end
fig, ax = plt.subplots()
ax.loglog(r_range, density_profile)
ax.loglog(rr, n0s)




streamlines = Qwind.interpolate_integrators(integrators, max_times=times);
mdot_sl_before = Qwind.compute_wind_mdot(streamlines, model.bh.Rg)
total_mass = 0
for loser_id in keys(merging)
    loser = streamlines[loser_id];
    winner = streamlines[merging[loser_id]["winner_id"]];
    if !escaped(loser)
        continue
    end
    winner_time = merging[loser_id]["winner_time"];
    loser_time = merging[loser_id]["loser_time"];
    winner_idx = searchsortedfirst(winner.t, winner_time)-1
    loser_idx = searchsortedfirst(loser.t, loser_time)-1
    v_loser = sqrt(loser.vr[loser_idx]^2 + loser.vz[loser_idx]^2)
    d_loser = sqrt(loser.r[loser_idx]^2 + loser.z[loser_idx]^2)
    d_winner_vec = @. sqrt(winner.r[winner_idx:end]^2 + winner.z[winner_idx:end]^2)
    v_winner_vec = @. sqrt(winner.vr[winner_idx:end]^2 + winner.vz[winner_idx:end]^2)
    loser_mass = 2π * loser.n[loser_idx] * d_loser * loser.width[loser_idx] * v_loser 
    v_winner = sqrt(winner.vr[winner_idx]^2 + winner.vz[winner_idx]^2)
    d_winner = sqrt(winner.r[winner_idx]^2 + winner.z[winner_idx]^2)
    winner_mass = 2π * winner.n[winner_idx] * d_winner * winner.width[winner_idx] * v_winner
    total_mass += loser_mass
    new_area = (loser.width[loser_idx] * d_loser + winner.width[winner_idx] * d_winner) 
    new_width = new_area / d_winner
    new_width = new_width * (d_winner_vec ./ d_winner)
    winner.width[winner_idx:end] .= new_width
    new_mass = winner_mass + loser_mass
    new_density = new_mass / (2π * d_winner * winner.width[winner_idx] * v_winner) # with the new width
    scaling_factor = d_winner^2 * v_winner * 1 ./ (d_winner_vec.^2 .* v_winner_vec)
    winner.n[winner_idx:end] = new_density * scaling_factor
end
escaping_sls = [streamline for streamline in streamlines if streamline.id ∉ collect(keys(merging))];
Qwind.compute_wind_mdot(escaping_sls, model.bh.Rg)
mdot_integ = Qwind.compute_wind_mdot(integrators, model.bh.Rg)
total_mass * C * M_P * 0.61 * model.bh.Rg^2
#println(mdot_integ - mdot_sl_before)

integrators = iterations_dict[1]["integrators"];
fig, ax = plt.subplots()
for integ in integrators
    #ax.plot(integ.p.data[:r], sqrt.(integ.p.data[:vr] .^2 + integ.p.data[:vz] .^2))
    ax.plot(integ.p.data[:r], integ.p.data[:z])
end
#ax.set_xlim(0,3000)
#ax.set_ylim(0,3000)

fig, ax = plt.subplots()
for integ in streamlines
    ax.plot(integ.r, integ.z)
end




getn0(model, 100)

compute_optical_depth(
    model.wind.density_grid,
    model.wind.grid_iterator,
    Qwind.ThomsonOpacity(),
    mu_electron = model.parameters.mu_electron,
    mu_nucleon = model.parameters.mu_nucleon,
    ri = 0,
    zi = 0,
    phii = 0,
    rf = 1000,
    zf = 1,
    phif = 0,
    source_luminosity = model.rad.xray_luminosity,
    Rg = model.bh.Rg,
    max_tau=Inf
)

#integs, sls = Qwind.run_integrators(model, iterations_dict; it_num=1, parallel = true);

#run_iteration!(model, iterations_dict, parallel = true, it_num=1);


@everywhere function compute_xi(rad, r, z; include_scattering = true)
    n = get_density(rad.wi.density_grid, r, z)
    tau_x = compute_tau_xray(model.wind.density_grid, model.wind.grid_iterator, model.rad, model.parameters, ri = 0, phii = 0, zi = 0, rf = r, zf = z, phif = 0)
    ret = compute_ionization_parameter(
        model.rad,
        model.wind,
        model.parameters,
        r = r,
        z = z,
        vr = 0,
        vz = 0,
        number_density = n,
        tau_x = tau_x,
    )
    return ret
end

@everywhere rr = collect(10 .^ range(log10(6), log10(1000), length = 50))
@everywhere zz = 10 .^ range(log10(1e-3), log10(1000), length = 51)
@everywhere rets = zeros(length(rr), length(zz))
@showprogress for (i, r) in enumerate(rr)
    rets[i, :] = pmap(
        z -> compute_xi(model.rad, r, z, include_scattering = true),
        zz,
        batch_size = 10,
    )
end

fig, ax = plt.subplots()
cm = ax.pcolormesh(rr, zz, rets', norm = LogNorm(), shading = "auto")
cbar = plt.colorbar(cm, ax = ax)
ax.set_xlabel("R [Rg]")
ax.set_ylabel("R [Rg]")
cbar.set_label("Ionization parameter", rotation = 270)




r = 1000
z = 1
phi = 0
n = 1e9
tau_x = compute_tau_xray(model.rad, ri = 0.0, phii = 0, zi = 0, rf = r, phif = phi, zf = z)
compute_ionization_parameter(model.rad, r = r, z = z, vr = 0, vz = 0, n = n, tau_x = tau_x)


rr = range(6, 1500, length = 500);
zz = 10 .^ range(-6, log10(1500), length = 500);
ret = zeros(length(rr), length(zz));
for (i, r) in enumerate(rr)
    for (j, z) in enumerate(zz)
        ret[i, j] = compute_tau_xray(
            model.rad,
            ri = 0.0,
            phii = 0.0,
            zi = model.rad.z_xray,
            rf = r,
            phif = 0.0,
            zf = z,
        )
    end
end

fig, ax = plt.subplots()
cm = ax.pcolormesh(rr, zz, ret', norm = LogNorm(vmin = 1e-2, vmax = 1e2))
plt.colorbar(cm, ax = ax)



scattered_grid = ScatteredLuminosityGrid()

rr = range(6, 1500, length = 500);
zz = 10 .^ range(-6, log10(1500), length = 500);
real_grid = zeros(length(rr), length(zz));
int_grid = zeros(length(rr), length(zz));
for (i, r) in enumerate(rr)
    for (j, z) in enumerate(zz)
        taux = compute_tau_xray(
            model.rad,
            ri = 0.0,
            phii = 0.0,
            zi = model.rad.z_xray,
            rf = r,
            phif = 0.0,
            zf = z,
        )
        density = interpolate_density(model.rad.wi.density_grid, r, z)
        xi = compute_ionization_parameter(model.rad, r, z, 1e-4, 1e-4, density, taux)
        int_grid[i, j] =
            interpolate_ionization_parameter(model.rad.wi.ionization_grid, r, z)
        real_grid[i, j] = xi
    end
end


fig, ax = plt.subplots(1, 2)
ax[1].pcolormesh(rr, zz, real_grid', norm = LogNorm())
ax[2].pcolormesh(rr, zz, int_grid', norm = LogNorm())

fig, ax = plt.subplots()
cm = ax.pcolormesh(rr, zz, (real_grid ./ int_grid)', norm = LogNorm())
plt.colorbar(cm, ax = ax)

iterations_dict[1] = Dict()
integrators, streamlines =
    Qwind.run_integrators!(model, iterations_dict, it_num = 1, parallel = true);
wi = model.rad.wi;

#new_radiation = update_radiation(model.rad, streamlines)

hull = Hull(streamlines, rtol = 1e-2);

density_grid = update_density_grid(wi.density_grid, wi.update_grid_flag, streamlines, hull);
velocity_grid =
    update_velocity_grid(wi.velocity_grid, wi.update_grid_flag, streamlines, hull);


#Profile.clear()
@time ionization_grid = IonizationGrid(
    density_grid,
    Rg = model.rad.bh.Rg,
    xray_luminosity = model.rad.xray_luminosity,
    z_xray = model.rad.z_xray,
    parallel = false,
);
#pprof()

fig, ax = plt.subplots()
cm = ax.pcolormesh(
    ionization_grid.r_range,
    ionization_grid.z_range,
    ionization_grid.grid',
    norm = LogNorm(),
)
plt.colorbar(cm, ax = ax)
#
#fig, ax = plt.subplots()
#cm = ax.pcolormesh(
#    density_grid.r_range,
#    density_grid.z_range,
#    density_grid.grid',
#    norm = LogNorm(),
#)
#plt.colorbar(cm, ax = ax)

ion0 = Qwind.compute_initial_ionization_grid(
    density_grid,
    model.rad.xray_luminosity,
    model.bh.Rg,
    model.rad.z_xray,
);
absorbed_from_center = compute_luminosity_absorbed_grid(
    density_grid,
    ion0,
    Rg = model.bh.Rg,
    source_position = [0, 0, model.rad.z_xray],
    source_luminosity = model.rad.xray_luminosity,
    mu_electron = 1.17,
    mu_nucleon = 0.61,
);


nodes, weights = gausslegendre(10);
i = 25;
j = 26;
#Profile.clear()
@time Qwind.compute_scattered_flux_in_cell(
    density_grid,
    ion0;
    scattered_luminosity_per_cell = absorbed_from_center,
    cell = Rectangle(
        density_grid.r_range[i],
        density_grid.r_range[i + 1],
        density_grid.z_range[j],
        density_grid.z_range[j + 1],
    ),
    cell_density = 1e8,
    absorption_opacity = BoostOpacity(),
    mu_electron = 1.17,
    mu_nucleon = 0.61,
    Rg = model.bh.Rg,
    nodes = nodes,
    weights = weights,
)
#pprof()

times = zeros(length(density_grid.r_range) - 1, length(density_grid.z_range) - 1);
for i = 1:(length(density_grid.r_range) - 1)
    for j = 1:(length(density_grid.z_range) - 1)
        r_source = (density_grid.r_range[i + 1] + density_grid.r_range[i]) / 2
        z_source = (density_grid.z_range[j + 1] + density_grid.z_range[j]) / 2
        time = @elapsed compute_optical_depth(
            density_grid.iterator,
            density_grid,
            ion0,
            BoostOpacity(),
            ri = 0,
            phii = π / 2,
            zi = 0,
            rf = density_grid.r_range[i],
            phif = 0,
            zf = density_grid.z_range[j],
            Rg = model.bh.Rg,
            mu_electron = 1.17,
            mu_nucleon = 0.61,
            max_tau = 30,
        )
        times[i, j] = time
    end
end

fig, ax = plt.subplots()
cm = ax.pcolormesh(
    density_grid.r_range[1:(end - 1)],
    density_grid.z_range[1:(end - 1)],
    times',
    norm = LogNorm(),
)
plt.colorbar(cm, ax = ax)
