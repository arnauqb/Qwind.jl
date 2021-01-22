using NearestNeighbors, Distances, Sundials
export create_wind_kdtree,
    get_density,
    WindKDTree,
    get_closest_points,
    get_closest_point,
    get_density_points,
    is_point_outside_wind

struct WindKDTree
    r::Vector{Float64}
    z::Vector{Float64}
    zmax::Vector{Float64}
    z0::Vector{Float64}
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
    z0_dense = Float64[]
    line_width_dense = Float64[]
    density_dense = Float64[]
    @info "Getting dense solution from integrators..."
    for integrator in integrators
        integration_time = unique(integrator.sol.t)
        tmin = integration_time[2]
        tmax = integration_time[end-1]
        if tmax <= tmin
            continue
        end
        t_range =
            10 .^ range(
                log10(tmin),
                log10(tmax),
                length = n_timesteps,
            )
        i = 1
        while (t_range[end] > integrator.sol.t[end]) && length(t_range) > 0
            t_range = t_range[1:end-i]
            i += 1
        end
        if length(t_range) == 0
            continue
        end
        dense_solution = integrator.sol(t_range)
        r_dense = vcat(r_dense, dense_solution[1, :])
        z_dense = vcat(z_dense, dense_solution[2, :])
        vr_dense = vcat(vr_dense, dense_solution[3, :])
        vz_dense = vcat(vz_dense, dense_solution[4, :])
        line_width_dense =
            vcat(line_width_dense, dense_solution[1, :] .* integrator.p.lwnorm)
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
        zmax_dense = vcat(
            zmax_dense,
            maximum(dense_solution[2, :]) .* ones(length(dense_solution[2, :])),
        )
        z0_dense = vcat(
            z0_dense,
            dense_solution[2, 1] .* ones(length(dense_solution[2, :])),
        )
    end
    @info "Done"
    return r_dense, z_dense, zmax_dense, z0_dense, line_width_dense, density_dense
end

function create_wind_kdtree(r, z, zmax, z0, width, density)
    points = hcat(r, z)'
    points = convert(Array{Float64,2}, points)
    kdtree = KDTree(points)
    wind_kdtree = WindKDTree(r, z, zmax, z0, width, density, kdtree)
    return wind_kdtree
end

function create_wind_kdtree(integrators::Array{Sundials.IDAIntegrator}, n_timesteps = 10000)
    r, z, zmax, z0, width, density=
        get_dense_solution_from_integrators(integrators, n_timesteps)
    return create_wind_kdtree(r, z, zmax, z0, width, density)
end

function get_closest_point(windkdtree::WindKDTree, r, z)
    idcs, distance = nn(windkdtree.tree, [r, z])
    rc = windkdtree.r[idcs]
    zc = windkdtree.z[idcs]
    zmaxc = windkdtree.zmax[idcs]
    z0c = windkdtree.z0[idcs]
    widthc = windkdtree.width[idcs]
    nc = windkdtree.n[idcs]
    return rc, zc, zmaxc, z0c, widthc, nc, distance
end

function get_closest_points(windkdtree::WindKDTree, points)
    idcs, distances = nn(windkdtree.tree, points)
    rc = windkdtree.r[idcs]
    zc = windkdtree.z[idcs]
    zmaxc = windkdtree.zmax[idcs]
    z0c = windkdtree.z0[idcs]
    widthc = windkdtree.width[idcs]
    nc = windkdtree.n[idcs]
    return rc, zc, zmaxc, z0c, widthc, nc, distances
end


function is_point_outside_wind(z, zmax, z0, distance, width)
    if (z > zmax) || (distance > width) || (z < z0)
        return true
    else
        return false
    end
end

"""
"""
function get_density(windkdtree::WindKDTree, r, z)
    rc, zc, zmax, z0, width, density, distance = get_closest_point(windkdtree, r, z)
    if is_point_outside_wind(z, zmax, z0, distance, width)
        return 1e2
    else
        return density
    end
end


function get_density_points(windkdtree::WindKDTree, points)
    rc, zc, zmax, z0, width, density, distances = get_closest_points(windkdtree, points)
    z = points[2, :]
    mask = is_point_outside_wind.(z, zmax, z0, distances, width) #(z .> zmax) || (distances .> width)
    density[mask] .= 1e2
    density
end
