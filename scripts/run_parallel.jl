using Distributed
@everywhere using DrWatson
@everywhere @quickactivate "Qwind"
#@everywhere using Qwind
using Qwind
using YAML
#include("scripts/plotting.jl")
#matplotlib.rcParams["figure.dpi"] = 300

config_path = "configs/multiple_models.yaml"
config = YAML.load_file(config_path, dicttype = Dict{Symbol,Any})
try
    mv(config[:integrator][:save_path], "backup", force = true)
catch
end
model = Model(config_path);

iterations_dict = Dict();
run!(model, iterations_dict)


do_iteration!(model, iterations_dict, it_num=1);

do_iteration!(model, iterations_dict, it_num=2);


fig, ax =plt.subplots()
for integ in iterations_dict[2]["integrators"]
    ax.plot(integ.p.data[:r], integ.p.data[:z])
end


using Dierckx, PyPlot, PyCall, ConcaveHull, DataFrames

LogNorm = matplotlib.colors.LogNorm

function reduce_line(integ)
    r, z, zmax, z0, width, n = get_dense_solution_from_integrator(integ, 10000);
    rs = [r[1]]
    zs = [z[1]]
    ns = [n[1]]
    for (rp, zp, np) in zip(r, z, n)
        # discard if they are very similar points
        #if isapprox([rp, zp], [rs[end], zs[end]], rtol = 0.55)
        #    continue
        #end
        #z_ratio = abs(zp / zs[end])
        if (zp > zs[end]) #&& (z_ratio < 0.95 || z_ratio > 1.05)
            push!(rs, rp)
            push!(zs, zp)
            push!(ns, np)
        end
    end
    return rs, zs, ns
end

integrators = iterations_dict[2]["integrators"];
rs = []
zs = []
ns = []
for integ in integrators[1:end]
    r, z, n = reduce_line(integ)
    rs = vcat(rs, r)
    zs = vcat(zs, z)
    ns = vcat(ns, n)
end
ns = log10.(ns);
rs_training = rs[1:1:end];#200000];
zs_training = zs[1:1:end];#200000];
ns_training = ns[1:1:end];#200000];
#rs_training = rs[200000:1:end];
#zs_training = zs[200000:1:end];
#ns_training = ns[200000:1:end];

df = DataFrame(:r=>rs, :z=>zs, :n=>ns)
dfr = round.(df, sigdigits=2)
idcs = nonunique(dfr[!,[:r, :z]])
dff = df[.!idcs, :]

#points = hcat(rs, zs, ns);
#points = round.(points, sigdigits=1);
#
#points_idx = unique(z -> points[z], 1:length(points), dims=1);

#rs_training = dff.r;
#zs_training = dff.z;
#ns_training = dff.n;
spl = Spline2D(rs_training, zs_training, ns_training, kx=1, ky=1, s) #2*length(rs_training));

ns_interp = evaluate(spl, rs, zs);
fig, ax = plt.subplots()
#ax.plot(ns, label="real")
#ax.plot(ns_interp, label="interp")
ax.plot(ns_interp ./ ns, label="ratio interp / real")
ax.legend()


#r_range = 10 .^ range(log10(56), log10(61), length=500);
#z_range = 10 .^ range(log10(0.045), log10(0.06), length=500);
r_range = 10 .^ range(log10(61), log10(200), length=250);
z_range = 10 .^ range(log10(1e-6), log10(20), length=250);
den_grid = zeros((length(r_range), length(z_range)));
for (i,r) in enumerate(r_range)
    for (j,z) in enumerate(z_range)
        den_grid[i,j] = 10 .^ evaluate(spl, r, z)
    end
end
den_grid_old = zeros((length(r_range), length(z_range)));
rt = iterations_dict[2]["radiative_transfer"];
for (i,r) in enumerate(r_range)
    for (j,z) in enumerate(z_range)
        den_grid_old[i,j] = get_density(rt, r, z);
    end
end
# scipy
scipy_interpolate = pyimport("scipy.interpolate")
r_range_grid = r_range .* ones(length(z_range))';
z_range_grid = z_range' .* ones(length(r_range));
den_scipy = 10 .^ scipy_interpolate.griddata((rs, zs), ns, (r_range_grid, z_range_grid), method="linear", fill_value = 2);
#scipy_interpolator = scipy_interpolate.LinearNDInterpolator(collect(zip(rs, zs)), ns, fill_value = 2);
#den_scipy = 10 .^ scipy_interpolator(r_range_grid, z_range_grid);

#den_scipy = 10 .^ den_scipy;

function get_value(r, z, r_range, z_range, den_grid)
    r_idx = searchsorted_nearest(r_range, r)
    z_idx = searchsorted_nearest(z_range, z)
    return den_grid[r_idx, z_idx]
end

fig, ax = plt.subplots()
cm = ax.pcolormesh(r_range, z_range, den_scipy', norm=LogNorm())
plt.colorbar(cm, ax=ax)
for integ in integrators
    ax.plot(integ.p.data[:r], integ.p.data[:z], color="white")
end
ax.set_xlim(r_range[1], r_range[end])
ax.set_ylim(z_range[1], z_range[end])
ax.set_xscale("log")
ax.set_yscale("log")


get_value(51.64, 1.69e-6, r_range, z_range, den_grid)

intt = integrators[3]
ri, zi, _, _, _, ni = get_dense_solution_from_integrator(intt, 10000);

fig, ax = plt.subplots()
cm = ax.pcolormesh(r_range, z_range, (den_scipy ./ den_grid)', norm=LogNorm(1e-2, 1e2), cmap="coolwarm")
#cm = ax.pcolormesh(r_range, z_range, den_scipy', norm=LogNorm())
for line in integrators
    ax.plot(line.p.data[:r], line.p.data[:z], color="white")
end
#ax.scatter(ri, zi, color="white", s=0.05)
plt.colorbar(cm, ax=ax)
ax.set_xlim(r_range[1], r_range[end])
ax.set_ylim(z_range[1], z_range[end])

fig, ax = plt.subplots()
semilogy(intt.p.data[:z], intt.p.data[:n])
ax.set_xlim(0.04, 0.07)

fig, ax = plt.subplots()
cm = ax.pcolormesh(r_range, z_range, den_grid_old', norm=LogNorm())
ax.scatter(ri, zi, color="white", s=0.05)
plt.colorbar(cm, ax=ax)
ax.set_xlim(r_range[1], r_range[end])
ax.set_ylim(z_range[1], z_range[end])


points = Vector{Float64}[]
integrators = iterations_dict[1]["integrators"];
#rmin = Inf
#rmax = 0
for integ in integrators
    #r, z, zmax, z0, width, n = get_dense_solution_from_integrator(integ, 10000);
    r = integ.p.data[:r]
    z = integ.p.data[:z]
    mxval, mxindx = findmax(z)
    #rmin = min(rmin, minimum(r))
    #rmax = max(rmax, maximum(r))
    push!(points, [r[mxindx], z[mxindx]])
    push!(points, [r[1], z[1]])
end
#pushfirst!(points, [rmin, 0.0])
#push!(points, [rmax, 0.0])
hull = concave_hull(points)


points2 = Vector{Float64}[]
for integ in integrators
    r = integ.p.data[:r]
    z = integ.p.data[:z]
    z_max = 0
    for (rp, zp) in zip(r, z)
        if zp > z_max 
            push!(points2, [rp, zp])
            z_max = zp
        end
    end
end
hull = concave_hull(points2)

r_min = minimum([minimum(integ.p.data[:r]) for integ in integrators]);
r_max = maximum([maximum(integ.p.data[:r]) for integ in integrators]);
r_range_p = range(r_min, r_max, step=1);
z_range_p = zeros(length(r_range_p));
for integ in integrators
    r, z, zmax, z0, width, n = get_dense_solution_from_integrator(integ, 10000);
    #r = integ.p.data[:r]
    #z = integ.p.data[:z]
    for (rp, zp) in zip(r, z)
        r_idx = searchsorted_nearest(r_range_p, rp)
        z_range_p[r_idx] = max(z_range_p[r_idx], zp)
    end
end
points3 = [[integ.p.r0, 0] for integ in integrators];
for (rp, zp) in zip(r_range_p, z_range_p)
    zp == 0 && continue
    push!(points3, [rp, zp])
end
hull = concave_hull(points3)


pp = hcat(points2...)
scatter(pp[1,:], pp[2,:])
pp = hcat(hull.vertices...)
scatter(pp[1,:], pp[2,:], color="red")

rr = range(0, 1500, length=500)
zz = range(0, 100, length=500)
ggrid = zeros((length(rr), length(zz)))
#polygon = AG.createpolygon(points3)
for (i,r) in enumerate(rr)
    for (j,z) in enumerate(zz)
        if in_hull([r, z], hull) #AG.contains(polygon, AG.createpoint(r,z))
            ggrid[i,j] = 1
        else
            ggrid[i,j] = 0
        end
    end
end
fig, ax = plt.subplots()
ax.pcolormesh(rr, zz, ggrid')
for integ in integrators
    ax.plot(integ.p.data[:r], integ.p.data[:z], color="white")
end
hullp = hcat(hull.vertices...)
#ax.scatter(hullp[1,:], hullp[2,:], color="red")
pp = hcat(points3...)
#ax.scatter(pp[1,:], pp[2,:], color = "blue", s=5)
ax.set_xlim(rr[1], rr[end])
ax.set_ylim(zz[1], zz[end])


using PyPlot
LogNorm = matplotlib.colors.LogNorm
rr = range(0, 1000, length=500)
zz = range(0, 25, length=500)
ggrid = zeros((length(rr), length(zz)))
for (i, r) in enumerate(rr)
    for (j, z) in enumerate(zz)
        ggrid[i,j] = get_density(model.rt, r, z)
    end
end
fig, ax = plt.subplots()
cm = ax.pcolormesh(rr, zz, ggrid', norm=LogNorm(1e7, 5e9))
plt.colorbar(cm, ax=ax)

fig, ax = plt.subplots()
rtg = model.rt.density_interpolator.grid
cm = ax.pcolormesh(rtg.r_range, rtg.z_range, rtg.grid', norm=LogNorm(1e7, 5e9))
plt.colorbar(cm, ax=ax)
