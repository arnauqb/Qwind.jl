using DrWatson
@quickactivate "Qwind"
using PyPlot
LogNorm = matplotlib.colors.LogNorm
using Qwind
using Profile
using PProf

model = Model("paper2/plot_nn.yaml");
iterations_dict = Dict();

do_iteration!(model, iterations_dict, it_num=1);

integrators = iterations_dict[1]["integrators"]
plot_streamlines(integrators)

lines_kdtrees = Qwind.create_lines_kdtrees(iterations_dict[1]["integrators"], n_timesteps=10000);

fig, ax = plt.subplots()
for line in lines_kdtrees
    ax.plot(line.r, line.z)
end

nr = 150;
nz = 151;
#r_range, z_range = Qwind.get_spatial_grid(lines_kdtrees, nr, nz);
r_range = range(10, 20, length=nr)
z_range = range(0.001, 0.2, length=nz)
#r_range = range(31.5, 33.5, length=nr)
#z_range = range(0.1, 0.2, length=nz)
function get_density_grid(lines_kdtrees, r_range, z_range)
    density_grid = zeros(Float64, nr, nz);
    for (i, r) in enumerate(r_range)
        for (j, z) in enumerate(z_range)
            density_grid[i, j] = get_density(lines_kdtrees, r, z);
        end
    end
    return density_grid
end
@time density_grid = get_density_grid(lines_kdtrees, r_range, z_range);

fig, ax= plt.subplots();
cm = ax.pcolormesh(r_range, z_range, density_grid', norm=LogNorm());
plt.colorbar(cm, ax=ax);
for line in lines_kdtrees
    ax.scatter(line.r, line.z, color="white", s=0.1)
end
ax.set_xlim(r_range[1], r_range[end])
ax.set_ylim(z_range[1], z_range[end])
