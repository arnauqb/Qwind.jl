using Distributed
@everywhere using DrWatson
@everywhere @quickactivate "Qwind"
@everywhere using Qwind
using PyPlot
LogNorm = matplotlib.colors.LogNorm
Normalize = matplotlib.colors.Normalize
using Profile, PProf

function read_table(file)
    df = CSV.read(file, DataFrame, header=false);
    return Matrix(df)
end

ry_values = "/home/arnau/Downloads/GX13p1_density_ps_05Ledd/np_1000.txt"
