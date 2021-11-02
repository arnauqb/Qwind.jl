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
Rg=bh.Rg
mu_electron=1.17
lumin = compute_bolometric_luminosity(bh)
dgrid = iterations_dict[3]["rad"].wi.density_grid;
#dgrid = DensityGrid("./debug/results.hdf5", 2);


#nr = 50
#nz = 51
#dgrid = DensityGrid(
#    range(0, 1000, length = nr),
#    range(0, 1000, length = nz),
#    1e10 .* ones(nr, nz),
#);

absorbed_from_center = compute_luminosity_absorbed_grid(
    dgrid,
    Rg = bh.Rg,
    source_position = [0, 0, 0],
    source_luminosity = lumin,
);

fig, ax = plt.subplots()
cm = ax.pcolormesh(
    dgrid.r_range[1:(end - 1)],
    dgrid.z_range[1:(end - 1)],
    absorbed_from_center',
    norm = LogNorm(vmin=1e32, vmax=1e34),
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



ret = compute_total_flux_grid(dgrid; Rg=Rg, mu_electron = mu_electron, source_luminosity=lumin)
