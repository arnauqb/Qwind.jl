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
model1 = Model(config_path);
iterations_dict1 = Dict();
do_iteration!(model1, iterations_dict1, it_num=1);

do_iteration!(model1, iterations_dict1, it_num=2);

#run!(model1, iterations_dict1)

#do_iteration!(model1, iterations_dict1, it_num=2)


#run!(model1, iterations_dict1)
#run!(model, iterations_dict, start_it=6, n_iterations=2)

config_path2 = "configs/config_test2.yaml"
config2 = YAML.load_file(config_path2, dicttype = Dict{Symbol,Any})
try
    mv(config2[:integrator][:save_path], "backup", force = true)
catch
end
model2 = Model(config_path2);
iterations_dict2 = Dict();
run!(model2, iterations_dict2)

do_iteration!(model, iterations_dict, it_num=1)

#do_iteration!(model, iterations_dict, it_num=2)

fig, ax = plot_density_grid(iterations_dict2[4]["radiative_transfer"].density_interpolator.grid)

function make_pdf(iterations_dict, fname)
    pdfile = pdfpages.PdfPages("$(fname)_lines.pdf")
    for iteration in sort(collect(keys(iterations_dict)))
        println("itereation $iteration")
        if "integrators" ∉ keys(iterations_dict[iteration])
            continue
        end
        integrators = iterations_dict[iteration]["integrators"]
        fig, ax = plot_streamlines(integrators)
        ax.set_title("Iteration $iteration")
        pdfile.savefig(fig)
        plt.close(fig)
    end
    pdfile.close()
    pdfile = pdfpages.PdfPages("$(fname)_grids.pdf")
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
end


make_pdf(iterations_dict1, "new")
make_pdf(iterations_dict2, "old")


# rt ratios
it_num = 4
grid1 = iterations_dict1[it_num]["radiative_transfer"].density_interpolator.grid;
grid2 = iterations_dict2[it_num]["radiative_transfer"].density_interpolator.grid;
fig, ax = plt.subplots();
cm = ax.pcolormesh(grid1.r_range, grid1.z_range, (grid1.grid ./ grid2.grid)', norm=LogNorm())


