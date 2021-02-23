using Distributed
@everywhere using DrWatson
@everywhere @quickactivate "Qwind"
#@everywhere using Qwind
using Qwind
using YAML
#include("scripts/plotting.jl")
#matplotlib.rcParams["figure.dpi"] = 300

config_path = "configs/config_test.yaml"
config = YAML.load_file(config_path, dicttype = Dict{Symbol,Any})
try
    mv(config[:integrator][:save_path], "backup", force = true)
catch
end
model = Model(config_path);
iterations_dict = Dict();
run!(model, iterations_dict)

#do_iteration!(model, iterations_dict, it_num=1);



pdfp = matplotlib.backends.backend_pdf.PdfPages
pdfile = pdfp("multipage_pdf.pdf")
for it in 1:length(iterations_dict)
    try
        println(it)
        fig, ax = plt.subplots()
        plot_streamlines(iterations_dict[it]["integrators"], fig, ax)
        pdfile.savefig(fig)
        plt.close("all")
    catch
        continue
    end
end
pdfile.close()

rt = iterations_dict[19]["radiative_transfer"]
r_range = range(6, 1000, length=100)
z_range = range(6, 1000, length=100)
tauxg = zeros((length(r_range), length(z_range)))
for (i,r ) in enumerate(r_range)
    for (j,z ) in enumerate(z_range)
        tauxg[i,j] = compute_xray_tau(rt, r, z)
    end
end
fig, ax = plt.subplots()
cm = ax.pcolormesh(r_range, z_range, tauxg', norm=LogNorm(vmin=1e-3, vmax=1e1))
plt.colorbar(cm, ax=ax)
fig.savefig("asd.pdf")






fig, ax = plt.subplots()
plot_streamlines(iterations_dict[21]["integrators"], fig, ax)
fig.savefig("asd.pdf")


lkdt = Qwind.create_lines_kdtrees(iterations_dict[1]["integrators"]);

using BenchmarkTools

vig = VIGrid(lkdt, 500, 50);


fig, ax = plt.subplots()
gg = vig.grid
ax.pcolormesh(gg.r_range, gg.z_range, gg.grid', norm=LogNorm());


r_range, z_range = get_spatial_grid(lkdt, 1000, 50);

gg = pmap(z -> get_density.(Ref(lkdt), r_range, z), z_range);
gg = hcat(gg...)







#########################

using JLD2, Profile, PProf
@load "big_rt.jld2" rt

iterations_dict[1] = Dict();
model.rt = rt;
iterations_dict[1]["radiative_transfer"] = rt;


do_iteration!(model, iterations_dict, it_num=1);

lines_range, lines_widths = get_initial_radii_and_linewidths(model.ic, model.bh.Rg);
i = 100
integrator = Qwind.create_and_run_integrator(
    model,
    linewidth = lines_widths[i],
    r0 = lines_range[i],
    atol = model.config[:integrator][:atol],
    rtol = model.config[:integrator][:rtol],
    line_id = i,
)


#run!(model, iterations_dict1)
#do_iteration!(model, iterations_dict, it_num=1);
#do_iteration!(model, iterations_dict, it_num=2);

Profile.clear()
@profile rt = update_radiative_transfer(model.rt, iterations_dict[1]["integrators"]);
pprof()

rt = iterations_dict[2]["radiative_transfer"];

#using JLD2
#@save "big_rt.jld2" rt
