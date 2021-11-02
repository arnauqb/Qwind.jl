using Distributed
@everywhere using DrWatson
@everywhere @quickactivate "Qwind"
@everywhere using Qwind
using PyPlot
LogNorm = matplotlib.colors.LogNorm
Normalize = matplotlib.colors.Normalize
using Profile, PProf

model = Model("configs/debug.yaml");
iterations_dict = Dict();
run!(model, iterations_dict);



bh = BlackHole(1e8 * M_SUN, 0.5, 0.0);
Rg = bh.Rg
mu_electron = 1.17
lumin = compute_bolometric_luminosity(bh)
#dgrid = iterations_dict[3]["rad"].wi.density_grid;

using HDF5
rr = h5read("./rin_20.hdf5", "iteration_049/density_grid/r");
zz = h5read("./rin_20.hdf5", "iteration_049/density_grid/z");
grid = h5read("./rin_20.hdf5", "iteration_049/density_grid/grid");
dgrid = DensityGrid(rr, zz, grid);

rgrid = 10 .^ range(log10(6), log10(5e3), length = 50);
zgrid = 10 .^ range(-6, log10(5e3), length = 50);
dgrid_interp = zeros(length(rgrid), length(zgrid));
for (i, r) in enumerate(rgrid)
    for (j, z) in enumerate(zgrid)
        dgrid_interp[i, j] = interpolate_density(dgrid, r, z)
    end
end
dgrid_interp = DensityGrid(rgrid, zgrid, dgrid_interp);

fig, ax = plt.subplots()
ax.pcolormesh(
    dgrid_interp.r_range,
    dgrid_interp.z_range,
    dgrid_interp.grid',
    norm = LogNorm(vmin = 1e2, vmax = 1e8),
)
#ax.set_xlim(0,1500)
#ax.set_ylim(0,1500)

#dgrid = DensityGrid("./rin_10.hdf5");


#nr = 50
#nz = 51
#dgrid = DensityGrid(
#    range(0, 1000, length = nr),
#    range(0, 1000, length = nz),
#    1e10 .* ones(nr, nz),
#);

absorbed_from_center = compute_luminosity_absorbed_grid(
    dgrid_interp,
    Rg = bh.Rg,
    source_position = [0, 0, 0],
    source_luminosity = lumin,
);

fig, ax = plt.subplots()
cm = ax.pcolormesh(
    dgrid_interp.r_range[1:(end - 1)],
    dgrid_interp.z_range[1:(end - 1)],
    absorbed_from_center',
    norm = LogNorm(),
)
plt.colorbar(cm, ax = ax)


cell = Rectangle(dgrid.r_range[7], dgrid.r_range[8], dgrid.z_range[7], dgrid.z_range[8])
compute_total_flux_in_cell(
    dgrid,
    scattered_luminosity_per_cell = absorbed_from_center,
    cell = cell,
    cell_density = 1e10,
    Rg = bh.Rg,
    mu_electron = 1.17,
    source_position = [0, 0, 0],
    source_luminosity = lumin,
)

#Profile.clear()
compute_scattered_flux_in_cell(
    dgrid,
    scattered_luminosity_per_cell = absorbed_from_center,
    cell = cell,
    cell_density = 1e12,
    Rg = bh.Rg,
    mu_electron = 1.17,
)
#pprof()



ret = compute_total_flux_grid(
    dgrid_interp;
    Rg = Rg,
    mu_electron = mu_electron,
    source_luminosity = lumin,
)

ret_center = Qwind.compute_total_flux_from_center_grid(
    dgrid_interp;
    Rg = Rg,
    mu_electron = mu_electron,
    source_luminosity = lumin,
);


fig, ax = plt.subplots(1, 2)
rr = dgrid_interp.r_range[1:(end - 1)]
zz = dgrid_interp.z_range[1:(end - 1)]
base_xi = zeros(length(rr), length(zz));
new_xi = zeros(length(rr), length(zz));
n_plot = zeros(length(rr), length(zz));
for (i, r) in enumerate(rr)
    for (j, z) in enumerate(zz)
        tau = compute_tau_uv(
            dgrid_interp,
            ri = 0,
            phii = 0,
            zi = 0,
            rf = r,
            zf = z,
            phif = 0,
            Rg = Rg,
            mu_electron = mu_electron,
        )
        n = interpolate_density(dgrid_interp, r, z)
        d = sqrt(r^2 + z^2) * Rg
        base_xi[i, j] = ret_center[i, j] / n#lumin / d^2 * exp(-tau) / n
        new_xi[i, j] = ret[i, j] / n
        n_plot[i, j] = n
    end
end
ax[1].pcolormesh(rr, zz, base_xi', norm = LogNorm(1, 1e10))
cm = ax[2].pcolormesh(rr, zz, new_xi', norm = LogNorm(1, 1e10))
plt.colorbar(cm, ax = ax[2])
for axis in ax
    axis.set_xlim(0, 1000)
    axis.set_ylim(0, 1000)
end


fig, ax = plt.subplots()
cm = ax.pcolormesh(rr, zz, (new_xi ./ base_xi)', norm=LogNorm())
plt.colorbar(cm, ax = ax)
ax.set_xscale("log")
ax.set_yscale("log")
