using Distributed
@everywhere using DrWatson
@everywhere @quickactivate "Qwind"
using Qwind
using YAML
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




do_iteration!(model, iterations_dict, it_num=1)

do_iteration!(model, iterations_dict, it_num=2)

plot_density_grid(iterations_dict[2]["radiative_transfer"].density_interpolator.grid, xlim=(0,250), ylim=(0,1))

lkdt = Qwind.create_lines_kdtrees(iterations_dict[1]["integrators"]);

r_range = ggrid.r_range #range(0, 250, length=100);
z_range = ggrid.z_range #range(0, 1, length=100);
lkdt_grid = zeros((length(r_range), length(z_range)));
for (i, r) in enumerate(r_range)
    for (j, z) in enumerate(z_range)
        lkdt_grid[i,j] = get_density(lkdt, r, z);
    end
end
fig, ax = plt.subplots();
cm = ax.pcolormesh(r_range, z_range, lkdt_grid', norm=LogNorm());
plt.colorbar(cm, ax=ax);



ggrid = iterations_dict[2]["radiative_transfer"].density_interpolator.grid;
lkgrid = zeros((length(ggrid.r_range), length(ggrid.z_range)));
for (i, r) in enumerate(ggrid.r_range)
    for (j, z) in enumerate(ggrid.z_range)
        lkgrid[i,j] = get_density(lkdt, r, z);
    end
end
ratio = ggrid.grid ./ lkgrid;
fig, ax = plt.subplots();
cm = ax.pcolormesh(ggrid.r_range, ggrid.z_range, ratio');
plt.colorbar(cm, ax=ax)










# Analysis

pdfile = pdfpages.PdfPages("lines.pdf")
for iteration in sort(collect(keys(iterations_dict)))
    integrators = iterations_dict[iteration]["integrators"]
    fig, ax = plot_streamlines(integrators)
    ax.set_title("Iteration $iteration")
    pdfile.savefig(fig)
    plt.close(fig)
end
pdfile.close()

pdfile = pdfpages.PdfPages("grids.pdf")
for iteration in sort(collect(keys(iterations_dict)))
    (iteration == 1) && continue
    rt = iterations_dict[iteration]["radiative_transfer"]
    grid = rt.density_interpolator.grid
    fig, ax = plot_density_grid(grid)
    ax.set_title("Iteration $iteration")
    pdfile.savefig(fig)
    plt.close(fig)
end
pdfile.close()
