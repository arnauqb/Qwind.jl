using NearestNeighbors
export create_wind_kdtree, get_density, WindKDTree

struct WindKDTree
    r::Vector{Float64}
    z::Vector{Float64}
    width::Vector{Float64}
    n::Vector{Float64}
    tree::KDTree
end

function create_wind_kdtree(
    solvers;
    n_time = 10000,
    tau_max_std = 0.01,
    cell_min_size = 100.0,
)
    r_dense = Float64[]
    z_dense = Float64[]
    n_dense = Float64[]
    v_r_dense = Float64[]
    v_z_dense = Float64[]
    line_width_dense = Float64[]
    density_dense = Float64[]
    for solver in solvers
        t_range =
            10 .^ range(
                max(log10(solver.sol.t[1]), 1e-6),
                log10(solver.sol.t[end-2]),
                length = n_time,
            )
        dense_solution = solver.sol(t_range)
        r_dense = vcat(r_dense, dense_solution[1, :])
        z_dense = vcat(z_dense, dense_solution[2, :])
        v_r_dense = vcat(v_r_dense, dense_solution[3, :])
        v_z_dense = vcat(v_z_dense, dense_solution[4, :])
        line_width_dense = vcat(
            line_width_dense,
            dense_solution[1, :] .* solver.p.line.width_norm,
        )
        v_t = @. sqrt(dense_solution[3, :]^2 + dense_solution[4, :]^2)
        density_dense = vcat(
            density_dense,
            compute_density.(
                dense_solution[1, :],
                dense_solution[2, :],
                v_t,
                Ref(solver.p.line),
            ),
        )

    kdtree = KDTree(points)
    wind_kdtree = WindKDTree(r_dense, z_dense, line_width_dense, density_dense, kdtree)
    return wind_kdtree
end

function getdensity(windkdtree::WindKDTree, r, z)
    idcs, distances = knn(windkdtree.tree, [r, z], 2)
    r1, r2 = windkdtree.r[idcs]
    z1, z2 = windkdtree.z[idcs]
    point = [r1, z1] + [r2-r1, z2-z1] / 2
    distance = evaluate(Euclidean(), point, [r,z])
    width = sum(windkdtree.width[idcs]) / 4 #2 of average and 2 of half width
    density = sum(windkdtree.n[idcs]) / 2
    if distance <= width #&& x < 0
        return density
    else
        return max(density * exp(-(distance - width)^2 / 2), 1e2)
    end
end
