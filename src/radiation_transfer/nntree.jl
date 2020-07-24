using NearestNeighbors, Distances
export create_wind_kdtree, get_density, WindKDTree

struct WindKDTree
    r::Vector{Float64}
    z::Vector{Float64}
    zmax::Vector{Float64}
    width::Vector{Float64}
    n::Vector{Float64}
    tree::KDTree
end

function create_wind_kdtree(
    integrators,
    n_time = 10000,
    tau_max_std = 0.001,
    cell_min_size = 0.001,
)
    r_dense = Float64[]
    z_dense = Float64[]
    n_dense = Float64[]
    vr_dense = Float64[]
    vz_dense = Float64[]
    zmax_dense = Float64[]
    line_width_dense = Float64[]
    density_dense = Float64[]
    for integrator in integrators
        t_range =
            10 .^ range(
                max(log10(integrator.sol.t[1]), 1e-6),
                log10(integrator.sol.t[end-2]),
                length = n_time,
            )
        dense_solution = integrator.sol(t_range)
        r_dense = vcat(r_dense, dense_solution[1, :])
        z_dense = vcat(z_dense, dense_solution[2, :])
        vr_dense = vcat(vr_dense, dense_solution[3, :])
        vz_dense = vcat(vz_dense, dense_solution[4, :])
        line_width_dense = vcat(
            line_width_dense,
            dense_solution[1, :] .* integrator.p.lwnorm,
        )
        v_t = @. sqrt(dense_solution[3, :]^2 + dense_solution[4, :]^2)
        density_dense = vcat(
            density_dense,
            compute_density.(
                dense_solution[1, :],
                dense_solution[2, :],
                dense_solution[3, :],
                dense_solution[4, :],
                Ref(integrator.p),
            ),
        )
        zmax_dense = vcat(zmax_dense, maximum(z_dense) .* ones(length(dense_solution[2,:])))
    end
    points = hcat(r_dense, z_dense)'
    kdtree = KDTree(points)
    wind_kdtree =
        WindKDTree(r_dense, z_dense, zmax_dense, line_width_dense, density_dense, kdtree)
    return wind_kdtree
end

function getdensity(windkdtree::WindKDTree, r, z)
    idcs, distances = knn(windkdtree.tree, [r, z], 2)
    r1, r2 = windkdtree.r[idcs]
    z1, z2 = windkdtree.z[idcs]
    maxz1, maxz2 = windkdtree.zmax[idcs]
    if z  > max(maxz1, maxz2)
        return 1e2
    end
    point = [r1, z1] + [r2 - r1, z2 - z1] / 2
    distance = evaluate(Euclidean(), point, [r, z])
    width = sum(windkdtree.width[idcs]) / 4 #2 of average and 2 of half width
    density = sum(windkdtree.n[idcs]) / 2
    if distance <= width #&& x < 0
        return density
    else
        return 1e2 #max(density * exp(-(distance - width)^2), 1e2)
    end
end
