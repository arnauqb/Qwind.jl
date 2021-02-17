using Distributed
@everywhere using DrWatson
@everywhere @quickactivate "Qwind"
@everywhere using Qwind
using YAML, HDF5, CSV, DataFrames
include("scripts/plotting.jl")
matplotlib.rcParams["figure.dpi"] = 300

mv(model.config[:integrator][:save_path], "backup",  force=true)
model = Model("results/resolution_tests/config1.yaml");
iterations_dict = Dict();
run!(model, iterations_dict)

xray_lumin = model.rad.xray_luminosity
Rg = model.bh.Rg

function compute_xray_grid(grid)
    r_range = grid["r"]
    z_range = grid["z"]
    vi_grid = VIGrid(r_range, z_range, grid["grid"])
    ret = zeros((length(r_range), length(z_range)))
    for (i,r) in enumerate(r_range)
        for (j, z) in enumerate(z_range)
            ret[i,j] = compute_xray_tau(vi_grid, 0.0, 0.0, r, z, xray_lumin, Rg)
            #ret[i,j] = compute_uv_tau(vi_grid, 0.0, 0.0, r, z, Rg)
        end
    end
    ret
end

function plot_streamlines(streamlines, fig, ax; color="white")
    line_ids = unique(streamlines.line_id)
    for lineid in line_ids
        mask = streamlines[!,:line_id] .== lineid
        toplot = streamlines[mask, :]
        ax.plot(toplot.r, toplot.z, color=color, linewidth=0.5)
    end
end

function plot_grid_lines(grid, fig, ax)
    for r in grid["r"]
        ax.axvline(r, color="black", linewidth=0.5)
    end
    for z in grid["z"]
        ax.axhline(z, color="black", linewidth=0.5)
    end
end


path = "results/resolution_tests/low_res"
lr_grid = h5read(path * "/results.hdf5", "iteration_002/density_grid");
lr_streamlines = CSV.read(path * "/iteration_002/streamlines.csv", DataFrame);
lr_vi_grid = VIGrid(lr_grid["r"], lr_grid["z"], lr_grid["grid"]);

hr_grid = h5read("results/auto_nlines/model_001/results.hdf5", "iteration_050/density_grid");
hr_vi_grid = VIGrid(hr_grid["r"], hr_grid["z"], hr_grid["grid"]);

fig, ax = plt.subplots();
tau_x = compute_xray_grid(lr_grid);
cm = ax.pcolormesh(lr_grid["r"], lr_grid["z"], tau_x', norm = LogNorm(vmin=1e-3));
plt.colorbar(cm, ax = ax);
plot_streamlines(lr_streamlines, fig, ax);
plot_grid_lines(lr_grid, fig, ax);
ax.set_xlim(50,54);
ax.set_ylim(0, 0.02);

fig, ax = plt.subplots()
tau_x = compute_xray_grid(hr_grid)
cm = ax.pcolormesh(, hr_z, tau_x', norm = Normalize(vmin=0.1))
plt.colorbar(cm, ax = ax)
#ax.set_xlim(2300, 2500)
ax.set_ylim(0.0, 1)

fig, ax = plt.subplots()
plot_streamlines(lr_streamlines, fig, ax, color="black")



r0 = 0
z0 = 0
r1 = 54
z1 = 0.02
m = (z1-z0) / (r1-r0)
n = z1 - m * r1
z(r) = m * r + n
r_range = range(r0, r1, length=10);
z_range = z.(r_range);
taux_range = compute_xray_tau.(Ref(lr_vi_grid), 0, 0, r_range, z_range, xray_lumin, Rg);
tauuv_range = compute_uv_tau.(Ref(lr_vi_grid), 0, 0, r_range, z_range, Rg);
plot(r_range, taux_range)
#plot(taux_range)

plot(r_range, tauuv_range)



density_grid = lr_vi_grid.interpolator(r_range, z_range);
fig, ax = plt.subplots()
cm = ax.pcolormesh(r_range, z_range, density_grid', norm = LogNorm())
plt.colorbar(cm, ax = ax)


gi = lr_vi_grid.grid.iterator
set_iterator!(gi, 0, 0, 54, 0.02)
points = copy(gi.intersection)
while !gi.finished
    next_intersection!(gi)
    points = hcat(points, gi.intersection)
end
points
