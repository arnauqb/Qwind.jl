using Qwind
using PyPlot
LogNorm = matplotlib.colors.LogNorm

it_num = 2
rt = iterations_dict[it_num]["radiative_transfer"]
sl = iterations_dict[it_num]["integrators"]

fig, ax = plt.subplots()
plot_streamlines(sl, fig, ax)

r_range = 10 .^ range(-1, 3, length=50)
z_range = 10 .^ range(-6, 3, length=50)

rfield = zeros((length(r_range), length(z_range), 2))
for (i, r) in enumerate(r_range)
    for (j, z) in enumerate(z_range)
        rfield[i, j, :] .= compute_disc_radiation_field(rt, r, z)
    end
end

force = compute_disc_radiation_field(rt, 6.250551925273973, 0.49417133613238345, atol=0, rtol=1e-3)

tau_uv = compute_uv_tau(rt, 28.96510171628526, 6.250551925273973, 0.49417133613238345)

fig, ax = plt.subplots()
cm = ax.pcolormesh(r_range, z_range, rfield[:,:,1]', norm=LogNorm())
plt.colorbar(cm, ax=ax)

fig, ax = plt.subplots()
cm = ax.pcolormesh(r_range, z_range, rfield[:,:,2]', norm=LogNorm(1e-10, 1e-5))
plt.colorbar(cm, ax=ax)

r_range = range(0, 200, length=250)
z_range = range(0, 100, length=250)
density_field = zeros((length(r_range), length(z_range)))
for (i, r) in enumerate(r_range)
    for (j, z) in enumerate(z_range)
        density_field[i, j] = get_density(rt.density_interpolator, r, z)
    end
end
fig, ax = plt.subplots()
cm = ax.pcolormesh(r_range, z_range, density_field', norm=LogNorm())
#plot_streamlines(sl, fig, ax)
plt.colorbar(cm, ax=ax)
ax.set_xlim(r_range[1], r_range[end])
ax.set_ylim(z_range[1], z_range[end])

fig, ax = plt.subplots()
cm = ax.pcolormesh(rt.density_interpolator.r_range, rt.density_interpolator.z_range, rt.density_interpolator.grid', norm=LogNorm())
#plot_streamlines(sl, fig, ax)
plt.colorbar(cm, ax=ax)
ax.set_xlim(0, 200)
ax.set_ylim(0, 1)
