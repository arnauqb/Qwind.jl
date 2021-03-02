using Distributed
@everywhere using DrWatson
@everywhere @quickactivate "Qwind"
using Qwind
using YAML, Profile, PProf, PyCall
include("scripts/plotting.jl")

config_path = "configs/config_test.yaml"
config = YAML.load_file(config_path, dicttype = Dict{Symbol,Any})
try
    mv(config[:integrator][:save_path], "backup", force = true)
catch
end
model1 = Model(config_path);
iterations_dict1 = Dict();
#do_iteration!(model1, iterations_dict1, it_num=1);
run!(model1, iterations_dict1)

integs = iterations_dict1[2]["integrators"];
r0 = [integ.p.r0 for integ in integs];
r, z, n = Qwind.reduce_integrators(integs, n_timesteps=100);

hull, points = Qwind.construct_wind_hull(r,z,r0);

vs = 10 .^ reduce(hcat, hull.vertices)
ps = 10 .^ reduce(hcat, points)

#fig, ax = plt.subplots()

fig, ax = QwindPlotting.plot_wind_hull(hull, nr=500, nz=500,rmin=49.999, rmax=50.01, zmax=10)
QwindPlotting.plot_streamlines(integs, fig, ax)
ax.set_xlim(49.999, 50.01)
ax.set_ylim(0, 10)
ax.scatter(vs[1,:], vs[2,:], color="blue", alpha=0.5)
#ax.scatter(ps[1,:], ps[2,:], color="red", alpha=0.5)
#ax.set_xlim(49.999, 50.003)
#ax.set_ylim(0, 1)
#ax.set_xlim(50.0, 50.03)
#ax.set_ylim(1e-6, 1)

r, z, n = Qwind.reduce_integrators(integs, n_timesteps=10000);

r2 = r[1:1:end];
z2 = z[1:1:end];
n2 = n[1:1:end];
log_n = log10.(n2);
r_range, z_range = get_spatial_grid(r2, z2, r0, "auto", 50);
r_range_grid = r_range .* ones(length(z_range))';
z_range_grid = z_range' .* ones(length(r_range));

scipy_interpolate = pyimport("scipy.interpolate");
#points = hcat(r, z);
#linear_int = scipy_interpolate.LinearNDInterpolator(points, log_n, fill_value=2)

#rbfi = scipy_interpolate.Rbf(r2, z2, log_n; function="linear");
using BasisFunctionExpansions

v = [log10.(r2) log10.(z2)];
rbf = MultiUniformRBFE(v, [10, 10], normalize=true)
bfa = BasisFunctionApproximation(log_n,v,rbf,0)

#vv = log10.([reduce(vcat, r_range_grid) reduce(vcat, z_range_grid)])
#yy = reshape(bfa(vv), length(z_range), length(r_range))

ret = zeros((length(r_range), length(z_range)));
for (i, r) in enumerate(r_range)
    for (j, z) in enumerate(z_range)
        if !Qwind.is_point_in_wind(hull, [r,z])
            ret[i,j] = 1e2
        else
            ret[i, j] = 10 .^ bfa(log10.([r z]))[1]
        end
    end
end

fig, ax = plt.subplots()
cm = ax.pcolormesh(r_range, z_range, ret', norm=LogNorm(vmin=1e5, vmax=1e10))
#QwindPlotting.plot_streamlines(integs, fig, ax, color="white", alpha=0.1)
ax.scatter(r[1:100:end], z[1:100:end], c=n[1:100:end], s=10, cmap="viridis", norm=LogNorm(vmin=1e5, vmax=1e10), edgecolor="black", linewidth=0.1)
plt.colorbar(cm, ax=ax)
ax.set_xlim(50, 50.05)
ax.set_ylim(0, 2)
#ax.set_xscale("log")


#yym = reshape(yy, (length(z_range), length(r_range)))

plot(r2, y)
plot!(r2, yy)

asd = yym'
asd = asd[abs.(asd) .!= Inf]
asd = reshape(asd, (length(z_range), length(r_range)))

heatmap(r_range, z_range, asd')


density_grid = 10 .^ rbfi(r_range_grid, z_range_grid)


#grid = Qwind.construct_interpolation_grid(r, z, n, r0, hull, nr="auto", nz=50);

#wi = WindInterpolator(integs);

#run!(model1, iterations_dict1)

#Profile.clear()
#@profile do_iteration!(model1, iterations_dict1, it_num=2);

xray_luminosity = model1.rad.xray_luminosity
Rg = model1.bh.Rg

# XRAY
it_num = 3
di = iterations_dict1[it_num]["radiative_transfer"].density_interpolator;
di2 = iterations_dict2[it_num]["radiative_transfer"].density_interpolator;
fig, ax = QwindPlotting.plot_xray_grid(di.grid, xray_luminosity, Rg, rmin=0, rmax=200, nr=250, nz=250, vmin=1e-1, vmax=1e4)
ax.set_yscale("log")
#ax.set_xscale("log")
ax.set_title("New")
fig, ax = QwindPlotting.plot_xray_grid(di2.grid, xray_luminosity, Rg, rmin=0, rmax=200, nr=250, nz=250, vmin=1e-1, vmax=1e4)
ax.set_title("Old")
ax.set_yscale("log")
#ax.set_xscale("log")




# DENSITY
it_num = 3
di = iterations_dict1[it_num]["radiative_transfer"].density_interpolator;
di2 = iterations_dict2[it_num]["radiative_transfer"].density_interpolator;
fig, ax = QwindPlotting.plot_density_grid(di.grid, xlim=(0,1e3), ylim=(1e-6,1e3))
ax.set_yscale("log")
#ax.set_xscale("log")
ax.set_title("New")
fig, ax = QwindPlotting.plot_density_grid(di2.grid, xlim=(0,1e3), ylim=(1e-6,1e3))
ax.set_title("Old")
ax.set_yscale("log")
#ax.set_xscale("log")

QwindPlotting.plot_streamlines(iterations_dict1[4]["integrators"])
QwindPlotting.plot_streamlines(iterations_dict2[4]["integrators"])


integ2 = iterations_dict2[3]["integrators"];
old_rt = iterations_dict2[4]["radiative_transfer"];
old_di = old_rt.density_interpolator;
new_rt = update_radiative_transfer(model1.rt, integ2);
new_di = new_rt.density_interpolator;

fig, ax = QwindPlotting.plot_xray_grid(new_di.grid, xray_luminosity, Rg, rmin=45, nr=500, nz=500, vmin=1e-2, vmax=1e2)
ax.set_yscale("log")
ax.set_xscale("log")
ax.set_title("New")
fig, ax = QwindPlotting.plot_xray_grid(old_di.grid, xray_luminosity, Rg, rmin=45, nr=500, nz=500, vmin=1e-2, vmax=1e2)
ax.set_title("Old")
ax.set_yscale("log")
ax.set_xscale("log")

it_num = 2
integ1 = iterations_dict1[it_num]["integrators"];
integ2 = iterations_dict2[it_num]["integrators"];
fig, ax = QwindPlotting.plot_streamlines(integ1)
ax.set_title("new")
fig, ax = QwindPlotting.plot_streamlines(integ2)
ax.set_title("old")
