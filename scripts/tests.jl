using Distributed
@everywhere using DrWatson
@everywhere @quickactivate "Qwind"
using YAML, Profile, PProf, PyCall
using Qwind
include("scripts/plotting.jl")

config_path = "configs/config_test.yaml"
config = YAML.load_file(config_path, dicttype = Dict{Symbol,Any})
try
    mv(config[:integrator][:save_path], "backup", force = true)
catch
end
model = Model(config_path);
iterations_dict = Dict();
run!(model, iterations_dict)

Profile.clear()
@profile do_iteration!(model, iterations_dict, it_num=2);

pprof()

vel_grid = model.rt.interpolator.velocity_grid;
den_grid = model.rt.interpolator.density_grid;
Rg = model.bh.Rg;

using PyPlot
LogNorm = matplotlib.colors.LogNorm
Normalize = matplotlib.colors.Normalize

fig, ax = plt.subplots()
for integ in iterations_dict[2]["integrators"]
    ax.plot(integ.p.data[:r], integ.p.data[:z], "o-")
end

function compute_momentum(r, z, vr, vz, n, Rg)
    θ = π / 2 - atan(z / r)
    sinθ = sin(θ)
    d = sqrt(r^2 + z^2)
    vt = sqrt(vr^2 + vz^2)
    momentum = vt^2.0 * n * 2 * π * d^2 * sinθ * Rg^2 * C^2 * M_P
    return momentum
end

thetas = π / 2 .- range(0, π/2, length=100);
r = 5000
R = r .* sin.(thetas);
z = r .* cos.(thetas);
n = interpolate_density.(Ref(den_grid), R, z);
vs = interpolate_velocity.(Ref(vel_grid), R, z);
vr = [v[1] for v in vs];
vz = [v[2] for v in vs];
momentums = compute_momentum.(R, z, vr, vz, n, Rg);
fig, ax = plt.subplots()
ax.plot(thetas, momentums)




r_range = collect(range(0, 5000, length = 100));
z_range = collect(range(0, 5000, length = 101));
r_range_grid = r_range .* ones(length(z_range))';
z_range_grid = z_range' .* ones(length(r_range));
vr_grid = zeros(length(r_range), length(z_range));
vz_grid = zeros(length(r_range), length(z_range));
density_grid_values = 1e2 .* ones(length(r_range), length(z_range));
momentum_grid = zeros(length(r_range), length(z_range));
for (i, r) in enumerate(r_range)
    for (j, z) in enumerate(z_range)
        vel = interpolate_velocity(vel_grid, r, z)
        vr_grid[i, j] = vel[1]
        vz_grid[i, j] = vel[2]
        n = interpolate_density(den_grid, r, z)
        density_grid_values[i, j] = n
        momentum_grid[i, j] = compute_momentum(r, z, vel[1], vel[2], n, Rg)
    end
end

fig, ax = plt.subplots()
cm = ax.pcolormesh(r_range, z_range, density_grid_values', norm = LogNorm())
ax.streamplot(
    r_range,
    z_range,
    vr_grid',
    vz_grid',
    density = 3,
    color="black"
    #color = momentum_grid',
    #norm = LogNorm(1e31, 5e34),
)
plt.colorbar(cm, ax = ax)
#ax.set_xlim(0,10000)
#ax.set_ylim(0,10000)

# calculate momentum at r = 5000 Rg

r = 5000
Rg = model.bh.Rg
thetas = range(0, π / 2, length = 20);
Rs = r .* sintheta;
zs = r .* cos.(thetas);
vs = interpolate_velocity.(Ref(vel_grid), Rs, zs);
vt = [sqrt(v[1]^2 + v[2]^2) for v in vs];
ns = interpolate_density.(Ref(den_grid), Rs, zs);
momentums = vt .^ 2.0 .* ns * 2 * π * r^2 .* sintheta * Rg^2 * C^2 * M_P;
fig, ax = plt.subplots()
ax.plot(90 .- (thetas * 180 / π), momentums)

QwindPlotting.plot_streamlines(iterations_dict[2]["integrators"])



# Check if integrators cross

using Interpolations
using Dierckx


integrators = iterations_dict[2]["integrators"];



integg1 = integrators[100];
integg2 = integrators[350];
fig, ax = plt.subplots()
ax.plot(integg1.p.data[:r], integg1.p.data[:z])
ax.plot(integg2.p.data[:r], integg2.p.data[:z])

using Roots


integ1 = DenseLine2(integg1);
integ2 = DenseLine2(integg2);

zcross = find_intersection(integ1, integ2)

fig, ax = plt.subplots()
ax.plot(integg1.p.data[:r], integg1.p.data[:z])
ax.plot(integg2.p.data[:r], integg2.p.data[:z])
ax.axhline(zcross, color = "black")

zmin = max(integ1.zmin, integ2.zmin)
zmax = min(integ1.zmax, integ2.zmax)
z_range = range(zmin, zmax, step = 1)

fig, ax = plt.subplots()
r1 = interpolate.(Ref(integ1), z_range)
r2 = interpolate.(Ref(integ2), z_range)
ax.plot(r1, z_range)
ax.plot(r2, z_range)

struct DenseLine10{T<:Vector{<:AbstractFloat}}
    r::T
    z::T
    vr::T
    vz::T
    n::T
    zmin::Float64
    zmax::Float64
    interpolator::Any
    function DenseLine10(integrator; n_timesteps = 1000)
        r, z, vr, vz, n = Qwind.reduce_line(integrator, n_timesteps = n_timesteps)
        interpolator = Interpolations.interpolate((z,), r, Gridded(Linear()))
        return new{typeof(r)}(r, z, vr, vz, n, minimum(z), maximum(z), interpolator)
    end
end

function interpolate(line::DenseLine10, z)
    if z < line.zmin
        return line.zmin
    elseif z > line.zmax
        return line.zmax
    else
        return line.interpolator(z)
    end
end

function find_intersection(integ1::DenseLine10, integ2::DenseLine10)
    zmin = max(integ1.zmin, integ2.zmin)
    zmax = min(integ1.zmax, integ2.zmax)
    (zmax < zmin) && return false
    z_range = range(zmin, zmax, length = 2000)
    f(z) = interpolate(integ1, z) - interpolate(integ2, z)
    previous_sign = sign(f(z_range[1]))
    i = 1
    for (i, z) in enumerate(z_range)
        current_sign = sign(f(z))
        if current_sign != previous_sign
            return z_range[i - 1]
        end
        previous_sign = current_sign
    end
    return Inf
end


function process_integrators(integrators)
    rs = Float64[]
    zs = Float64[]
    vrs = Float64[]
    vzs = Float64[]
    ns = Float64[]
    dense_lines = DenseLine10.(integrators, n_timesteps=10000)
    for (i, line) in enumerate(dense_lines)
        zint = Inf
        for (j, line2) in enumerate(dense_lines)
            if j <= i
                continue
            end
            zint = min(zint, find_intersection(line, line2))
        end
        z_idx = searchsorted_nearest(line.z, zint)
        rs = vcat(rs, line.r[1:z_idx])
        zs = vcat(zs, line.z[1:z_idx])
        vrs = vcat(vrs, line.vr[1:z_idx])
        vzs = vcat(vzs, line.vz[1:z_idx])
        ns = vcat(ns, line.n[1:z_idx])
    end
    return rs, zs, vrs, vzs, ns
end

struct Line
    r::Any
    z::Any
    vr::Any
    vz::Any
    n::Any
end

function process_integrators_individual(integrators)
    lines = Line[]
    dense_lines = DenseLine10.(integrators, n_timesteps=10000)
    for (i, line) in enumerate(dense_lines)
        zint = Inf
        for (j, line2) in enumerate(dense_lines)
            if j <= i
                continue
            end
            zint = min(zint, find_intersection(line, line2))
        end
        z_idx = searchsorted_nearest(line.z, zint)
        p_line = Line(
            line.r[1:z_idx],
            line.z[1:z_idx],
            line.vr[1:z_idx],
            line.vz[1:z_idx],
            line.n[1:z_idx],
        )
        push!(lines, p_line)
    end
    return lines
end

processed_integrators = process_integrators(integrators);

integs_indiv = process_integrators_individual(integrators);

fig, ax = plt.subplots()
for integ in integs_indiv
    ax.plot(integ.r, integ.z)
end

fig, ax = plt.subplots()
ax.scatter(processed_integrators[1][1:10:end], processed_integrators[2][1:10:end])

density_grid = Qwind.construct_density_grid(
    processed_integrators[1],
    processed_integrators[2],
    processed_integrators[5],
    r0s,
    hull,
    nr = "auto",
    nz = 500,
);

fig, ax = plt.subplots()
cm = ax.pcolormesh(
    density_grid.r_range,
    density_grid.z_range,
    density_grid.grid',
    norm = LogNorm(),
    shading = "auto",
)
plt.colorbar(cm, ax = ax)
QwindPlotting.plot_streamlines(
    integrators,
    fig = fig,
    ax = ax,
    color = "white",
    alpha = 0.1,
)
ax.set_ylim(1e-6, density_grid.z_range[end])

velocity_grid = Qwind.construct_velocity_grid(
    processed_integrators[1],
    processed_integrators[2],
    processed_integrators[3],
    processed_integrators[4],
    r0s,
    hull,
    nr = "auto",
    nz = 500,
    log = true,
);

r_range = collect(range(0, 5000, length = 50));
z_range = collect(range(0, 5000, length = 51));
r_range_grid = r_range .* ones(length(z_range))';
z_range_grid = z_range' .* ones(length(r_range));
vr_grid = zeros(length(r_range), length(z_range));
vz_grid = zeros(length(r_range), length(z_range));
density_grid_values = 1e2 .* ones(length(r_range), length(z_range));
momentum_grid = zeros(length(r_range), length(z_range));
for (i, r) in enumerate(r_range)
    for (j, z) in enumerate(z_range)
        vel = interpolate_velocity(velocity_grid, r, z)
        vr_grid[i, j] = vel[1]
        vz_grid[i, j] = vel[2]
        n = interpolate_density(density_grid, r, z)
        density_grid_values[i, j] = n
        momentum_grid[i, j] = compute_momentum(r, z, vel[1], vel[2], n, Rg)
    end
end
vt_grid = sqrt.(vr_grid .^2 + vz_grid .^2);

fig, ax = plt.subplots()
cm = ax.pcolormesh(
    r_range,
    z_range,
    density_grid_values',
    norm = LogNorm(vmax=1e6),
    shading = "auto",
)
ax.streamplot(
    r_range,
    z_range,
    vr_grid',
    vz_grid',
    density = 3,
    color = vt_grid',
    cmap = plt.get_cmap("spring"),
    norm=LogNorm(1e-2, 0.25),
    #color = momentum_grid',
    #norm = LogNorm(1e31, 5e34),
)
plt.colorbar(cm, ax = ax)
ax.set_ylim(10, 5000)
#ax.set_ylim(0,10000)


# compare grids
r_range = 10 .^ range(1, 4, length = 250);
z_range = 10 .^ range(-6, 4, length = 250);
r_range_grid = r_range .* ones(length(z_range))';
z_range_grid = z_range' .* ones(length(r_range));

original_grid = iterations_dict[3]["radiative_transfer"].interpolator.density_grid;

original_values = get_density.(Ref(original_grid), r_range_grid, z_range_grid);
new_values = get_density.(Ref(density_grid), r_range_grid, z_range_grid);

cmap = plt.get_cmap("RdBu")
fig, ax = plt.subplots()
cm = ax.pcolormesh(
    r_range,
    z_range,
    (new_values ./ original_values)',
    norm = LogNorm(vmin = 1e-1, vmax = 1e1),
    cmap = cmap,
)
plt.colorbar(cm, ax = ax)


hull_dense = Qwind.construct_wind_hull(processed_integrators[1], processed_integrators[2], r0s, sigdigits=10)


r_range = 10 .^ range(log10(1), log10(2e4), length=150);
z_range = 10 .^ range(log10(1e-6), log10(2e4), length=151);
old_hull_grid = zeros((length(r_range), length(z_range)));
new_hull_grid = zeros((length(r_range), length(z_range)));
for (i, r) in enumerate(r_range)
    for (j, z) in enumerate(z_range)
        old_hull_grid[i, j] = Qwind.is_point_in_wind(hull, r, z)
        new_hull_grid[i, j] = Qwind.is_point_in_wind(hull_dense, r, z)
    end
end

fig, ax = plt.subplots(1, 2, figsize=(10, 4))
ax[1].pcolormesh(r_range, z_range, old_hull_grid')
ax[2].pcolormesh(r_range, z_range, new_hull_grid')
QwindPlotting.plot_streamlines(integrators, fig=fig, ax=ax[1], color="white")
QwindPlotting.plot_streamlines(integrators, fig=fig, ax=ax[2], color="white")
for axis in ax
    axis.set_xlim(r_range[1], r_range[end])
    axis.set_ylim(z_range[1], z_range[end])
end



