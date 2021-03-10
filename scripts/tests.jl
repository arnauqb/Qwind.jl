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
model = Model(config_path);
iterations_dict = Dict();
run!(model, iterations_dict)


# compute momentum

using CSV, DataFrames, PyPlot

df = CSV.read("./runs/tests/iteration_002/streamlines.csv", DataFrame);

fig, ax = plt.subplots()
for line_id in unique(df.line_id)
    toplot = filter(row -> row.line_id == line_id, df)
    ax.plot(toplot.r, toplot.z, label=nothing)
end

total_momentum = 0.0
for line_id in unique(df.line_id)
    toplot = filter(row -> row.line_id == line_id, df)
    total_momentum += toplot.momentum[end]
end
