using NearestNeighbors, Distances, Sundials
export create_wind_kdtree, get_density, WindKDTree, get_closest_points

struct WindKDTree
    r::Vector{Float64}
    z::Vector{Float64}
    zmax::Vector{Float64}
    width::Vector{Float64}
    n::Vector{Float64}
    tree::KDTree
end

function get_dense_solution_from_integrators(integrators, n_timesteps = 10000)
    r_dense = Float64[]
    z_dense = Float64[]
    vr_dense = Float64[]
    vz_dense = Float64[]
    zmax_dense = Float64[]
    line_width_dense = Float64[]
    density_dense = Float64[]
    for integrator in integrators
        t_range =
            10 .^ range(
                max(log10(integrator.sol.t[1]), 1e-6),
                log10(integrator.sol.t[end - 2]),
                length = n_timesteps,
            )
        dense_solution = integrator.sol(t_range)
        r_dense = vcat(r_dense, dense_solution[1, :])
        z_dense = vcat(z_dense, dense_solution[2, :])
        vr_dense = vcat(vr_dense, dense_solution[3, :])
        vz_dense = vcat(vz_dense, dense_solution[4, :])
        line_width_dense =
            vcat(line_width_dense, dense_solution[1, :] .* integrator.p.lwnorm)
        v_t = @. sqrt(dense_solution[3, :] .^ 2 + dense_solution[4, :]^2)
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
        zmax_dense =
            vcat(zmax_dense, maximum(z_dense) .* ones(length(dense_solution[2, :])))
    end
    return r_dense, z_dense, zmax_dense, line_width_dense, density_dense
end

function create_wind_kdtree(r, z, zmax, width, density)
    points = hcat(r, z)'
    points = convert(Array{Float64,2}, points)
    kdtree = KDTree(points)
    wind_kdtree = WindKDTree(r, z, zmax, width, density, kdtree)
    return wind_kdtree
end

function create_wind_kdtree(integrators::Array{Sundials.IDAIntegrator}, n_timesteps = 10000)
    r, z, zmax, width, density =
        get_dense_solution_from_integrators(integrators, n_timesteps)
    return create_wind_kdtree(r, z, zmax, width, density)
end

function get_closest_points(windkdtree::WindKDTree, r, z, k = 1)
    idcs, distances = knn(windkdtree.tree, [r, z], k)
    rc = windkdtree.r[idcs]
    zc = windkdtree.z[idcs]
    zmaxc = windkdtree.zmax[idcs]
    widthc = windkdtree.width[idcs]
    nc = windkdtree.n[idcs]
    return rc, zc, zmaxc, widthc, nc
end


"""
Gets density at the point [r,z] using the windtree WindTree.
We first get the two closest points belonging to a streamline to [r,z].
Then we compute the point halfway between them. If the distance to the point
is smaller than the width of the line at those points, then it returns the mean
density between the two points. Otherwise, it returns an exponential decay.
"""
function get_density(windkdtree::WindKDTree, r, z)
    (r1, r2), (z1, z2), (zmax1, zmax2), (width1, width2), (n1, n2) =
        get_closest_points(windkdtree, r, z, 2)
    if z > max(zmax1, zmax2)
        return 1e2
    end
    point = [r1, z1] + [r2 - r1, z2 - z1] / 2
    distance = evaluate(Euclidean(), point, [r, z])
    width = (width1 + width2) / 4 # 2 of average and 2 of half width
    density = (n1 + n2) / 2
    if distance <= width #&& x < 0
        return density
    else
        return 1e2 #max(density * exp(-(distance - width)^2), 1e2)
    end
end
