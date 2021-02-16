using Distributed
@everywhere using DrWatson
@everywhere @quickactivate "Qwind"
@everywhere using Qwind
using YAML
include("scripts/plotting.jl")
matplotlib.rcParams["figure.dpi"] = 300

# model 1
model1 = Model("results/resolution_tests/config1.yaml");
iterations_dict1 = Dict();
run!(model1, iterations_dict1)
#do_iteration!(model1, iterations_dict1, it_num=1);

# model 2
model2 = Model("results/resolution_tests/config2.yaml");
iterations_dict2 = Dict();
run!(model2, iterations_dict2)
#do_iteration!(model2, iterations_dict2, it_num=1);

# second iteration
#do_iteration!(model1, iterations_dict1, it_num=2);
#do_iteration!(model2, iterations_dict2, it_num=2);
#

angles = []
momentum = []
for line in iterations_dict1[2]["integrators"]
    angle = atan(line.p.data[:z][end] / line.p.data[:r][end])
    push!(angles, angle / π * 180)
    vel = sqrt(line.p.data[:vr][end]^2 + line.p.data[:vz][end]^2)
    mom = Qwind.compute_integrator_mdot(line, model1.bh.Rg) * vel * C
    push!(momentum, mom)
end
angles2 = range(0, 90, step=5)
moms2 = zeros(length(angles2))
for (i, (angle, mom)) in enumerate(zip(angles, momentum))
    idx = searchsorted_nearest(angles2, angle)
    moms2[idx] += mom
end
fig, ax = plt.subplots(figsize=(10,5))
ax.plot(angles2, moms2, "o-")
ax.set_xlabel("Angle [deg]")
ax.set_ylabel("Momentum / time [g cm / s^2]")
ax.set_yscale("log")
fig.savefig("Momentum_angle.png", dpi=300, bbox_inches="tight")



vt(vr, vz) = sqrt.(vr .^ 2 + vz .^ 2)
d(r, z) = sqrt.(r .^ 2 + z .^ 2)
function compute_xray_grid(rt, r_range, z_range)
    ret = zeros((length(r_range), length(z_range)))
    for (i,r) in enumerate(r_range)
        for (j, z) in enumerate(z_range)
            ret[i,j] = compute_xray_tau(rt, r, z)
        end
    end
    ret
end
function compute_density_grid(rt, r_range, z_range)
    ret = zeros((length(r_range), length(z_range)))
    for (i,r) in enumerate(r_range)
        for (j, z) in enumerate(z_range)
            ret[i,j] = get_density(rt, r, z)
        end
    end
    ret
end
function plot_grid(r_range, z_range, ax)
    for r in r_range
        ax.axvline(r, color="black", alpha=0.5)
    end
    for z in z_range
        ax.axhline(z, color="black", alpha=0.5)
    end
end

it_num = 3
rt = iterations_dict1[it_num]["radiative_transfer"];
r_range = 10 .^ range(-2, 3, length=100);
z_range = 10 .^ range(-3, 3, length=100);
den_grid = compute_xray_grid(rt, r_range, z_range);
fig, ax = plt.subplots();
cm = ax.pcolormesh(r_range, z_range, den_grid', norm=LogNorm(vmin=5e-2, vmax=1e-1));
#plot_streamlines(iterations_dict1[it_num]["integrators"], fig, ax)
ax.set_xlim(r_range[1], r_range[end])
ax.set_ylim(z_range[1], z_range[end])
ax.set_xlabel("R [Rg]")
ax.set_ylabel("z [Rg]")
ax_cbar = plt.colorbar(cm, ax=ax)
ax_cbar.set_label("Tau X")
ax.set_title("Bad shielding")
fig.savefig("xray_grid2.png", dpi=300, bbox_inches="tight")







fig, ax = plt.subplots()
for line in iterations_dict1[2]["integrators"]
    ax.loglog(line.p.data[:z], vt(line.p.data[:vr], line.p.data[:vz]), color="C0")
end
for line in iterations_dict2[2]["integrators"]
    ax.loglog(line.p.data[:z], vt(line.p.data[:vr], line.p.data[:vz]), color="C1")
end

line_id = 2
ll1 = iterations_dict1[2]["integrators"][line_id]
ll2 = iterations_dict2[2]["integrators"][line_id]
loglog(ll1.p.data[:r], ll1.p.data[:z])
loglog(ll2.p.data[:r], ll2.p.data[:z])

function max_v(integrator)
    vt = @. sqrt(integrator.p.data[:vr]^2 + integrator.p.data[:vz]^2)
    return maximum(vt)
end

it_num = 2
integ1 = iterations_dict1[it_num]["integrators"];
integ2 = iterations_dict2[it_num]["integrators"];
vs1 = max_v.(integ1);
vs2 = max_v.(integ2);
fig, ax = plt.subplots()
ax.plot(vs1, "o-", label="low res")
ax.plot(vs2, "o-", label="high res")
ax.legend()



plot_grid(rt, ax) = plot_grid(rt.density_interpolator.grid.r_range, rt.density_interpolator.grid.z_range, ax)

it_num = 2
rt1 = iterations_dict1[it_num]["radiative_transfer"];
rt2 = iterations_dict2[it_num]["radiative_transfer"];
integ1 = iterations_dict1[it_num-1]["integrators"]
integ2 = iterations_dict2[it_num-1]["integrators"]
r_range = range(50, 57, length=150);
z_range = range(1e-6, 20, length=150);
xray_grid1 = compute_xray_grid(rt1, r_range, z_range);
xray_grid2 = compute_xray_grid(rt2, r_range, z_range);
density_grid1 = compute_density_grid(rt1, r_range, z_range);
density_grid2 = compute_density_grid(rt2, r_range, z_range);

# xray
fig, ax= plt.subplots(1, 2)
cm = ax[1].pcolormesh(r_range, z_range, xray_grid1', norm=LogNorm(1e-2,1e2));
cm = ax[2].pcolormesh(r_range, z_range, xray_grid2', norm=LogNorm(1e-2,1e2));
ax[1].set_title("low res")
ax[2].set_title("high res")
plt.colorbar(cm, ax=ax[2])
for axis in ax
    plot_streamlines(integ1, fig, axis, color="blue")
    plot_streamlines(integ2, fig, axis, color="red")
    axis.set_xlim(r_range[1], r_range[end])
    axis.set_ylim(z_range[1], z_range[end])
end
plot_grid(rt1, ax[1])
plot_grid(rt2, ax[2])

# density
fig, ax= plt.subplots(1, 2)
cm = ax[1].pcolormesh(r_range, z_range, density_grid1', norm=LogNorm(vmin=1e6));
cm = ax[2].pcolormesh(r_range, z_range, density_grid2', norm=LogNorm(vmin=1e6));
ax[1].set_title("low res")
ax[2].set_title("high res")
#plt.colorbar(cm, ax=ax[2])
for axis in ax
    plot_streamlines(integ1, fig, axis, color="blue")
    plot_streamlines(integ2, fig, axis, color="red")
    axis.set_xlim(r_range[1], r_range[end])
    axis.set_ylim(z_range[1], z_range[end])
end

fig, ax = plt.subplots()
ratio = xray_grid2 ./ xray_grid1
cm = ax.pcolormesh(r_range, z_range, ratio', norm=LinearNorm(0.8, 1.1))
plt.colorbar(cm, ax=ax)

lines_kdtrees = Qwind.create_lines_kdtrees(iterations_dict1[1]["integrators"]);
lkd_grid = zeros((length(r_range), length(z_range)));
for (i,r) in enumerate(r_range)
    for (j, z) in enumerate(z_range)
        lkd_grid[i,j] = get_density(lines_kdtrees, r, z);
    end
end
fig, ax = plt.subplots(1, 2)
ggrid = rt1.density_interpolator.grid;
itp = interpolate((ggrid.r_range, ggrid.z_range), ggrid.grid, Gridded(Linear()));
intep_grid1 = zeros((length(r_range), length(z_range)));
for (i,r) in enumerate(r_range)
    for (j, z) in enumerate(z_range)
        intep_grid1[i,j] = itp(r, z)
    end
end
ratio1 = density_grid1 ./ lkd_grid;
ratio2 = density_grid2 ./ lkd_grid;
cm = ax[1].pcolormesh(r_range, z_range, ratio1', norm=LinearNorm(0.95, 1.05), cmap="coolwarm")
cm = ax[2].pcolormesh(r_range, z_range, ratio2', norm=LinearNorm(0.95, 1.05), cmap="coolwarm")
plt.colorbar(cm, ax=ax)
#for r in rt1.density_interpolator.grid.r_range
#    ax[1].axvline(r, color="black", alpha=0.5)
#end
#for z in rt1.density_interpolator.grid.z_range
#    ax[1].axhline(z, color="black", alpha=0.5)
#end
#for r in rt2.density_interpolator.grid.r_range
#    ax[2].axvline(r, color="black", alpha=0.5)
#end
#for z in rt2.density_interpolator.grid.z_range
#    ax[2].axhline(z, color="black", alpha=0.5)
#end
for axis in ax
    axis.set_xlim(r_range[1], r_range[end])
    axis.set_ylim(z_range[1], z_range[end])
end


# simple x-ray test
function make_grid(nr, nz)
    rr = range(50, 1000, length=nr);
    zz = range(0, 100, length=nz);
    grid = 1e8 .* ones((length(rr), length(zz)));
    vi_grid1 = VIGrid(collect(rr), collect(zz), grid);
end

gg0 = make_grid(2,2);
gg1 = make_grid(10, 2);
gg2 = make_grid(100, 2);
gg3 = make_grid(10000, 2);

xray_lumin = model1.rad.xray_luminosity
Rg = model1.bh.Rg

Qwind.compute_xray_tau_cell(500 * Rg, 500*Rg, 1e8, 0, xray_lumin, Rg)

r = 600
Qwind.compute_xray_tau_cell((r-50) * Rg, 50*Rg, 1e8, 0, xray_lumin, Rg)
taux1 = compute_xray_tau(gg1, 0, 0, r, 0, xray_lumin, Rg)

taux2 = compute_xray_tau(gg2, 0, 0, r, 0, xray_lumin, Rg)

taux3 = compute_xray_tau(gg3, 0, 0, r, 0, xray_lumin, Rg)


rx_range = range(50, 1000, length=50);
values = Qwind.ionization_cell_xi_kernel.(xray_lumin, 1e8, 50*Rg, 0, rx_range.*Rg);
plot(rx_range, values)
axhline(0)

using Optimize 

d0 = find_zero(t->Qwind.ionization_cell_xi_kernel(xray_lumin, 1e8, 50*Rg, 0, t), (50Rg, 550*Rg))
println("correct $(d0/Rg)")

r_range = range(10, 1000, length=50);
rd_range = range(10, 1000, length=50);
z_range = range(1e-6, 20, length=50);
uvgrid = zeros((length(r_range), length(z_range)));
rt = iterations_dict2[2]["radiative_transfer"];

function compute_uv_grid(rt, r_range, z_range, rd_range)
    for (i,r) in enumerate(r_range)
        for (j, z) in enumerate(z_range)
            for (k, rd) in enumerate(rd_range)
                uvgrid[i,j] = compute_uv_tau(rt, rd, r, z);
            end
        end
    end
    return uvgrid
end

using BenchmarkTools, Profile, PProf

Profile.clear()
@profile compute_uv_grid(rt, r_range, z_range, rd_range);
pprof()
