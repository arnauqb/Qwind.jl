using CSV, DataFrames

df = CSV.read("./results_kernel_size_0/iteration_006.csv", DataFrame)

kdtree = create_wind_kdtree(df.r, df.z, df.zmax, df.z0, df.width, df.density, 1e2, 100000)
