using DrWatson
@quickactivate "Qwind"
using RegionTrees, DataFrames, CSV, YAML, Printf, DataFrames
using Qwind

kdtree_data = CSV.read("kdtree_data.csv", DataFrame)

bh = BlackHole(1e8, 0.5, 0.5)

wind_kdtree = create_wind_kdtree(
    kdtree_data.r,
    kdtree_data.z,
    kdtree_data.zmax,
    kdtree_data.z_0,
    kdtree_data.width,
    kdtree_data.n,
)

quadtree = create_and_refine_quadtree(
    wind_kdtree,
    Rg = bh.Rg,
    atol = 1e-2,
    rtol = 1e-2,
    cell_min_size = 1e-6,
)
