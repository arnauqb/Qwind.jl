import NearestNeighbors
using Distances, Sundials
export create_wind_kdtree,
    create_line_kdtree,
    create_lines_kdtree,
    get_density,
    get_density_closest,
    r_min,
    r_max,
    z_min,
    z_max,
    KDTree,
    LineKDTree,
    get_closest_points,
    get_closest_point,
    get_density_points,
    is_point_outside_wind


struct KDTree <: NNInterpolator
    r::Vector{Float64}
    z::Vector{Float64}
    zmax::Vector{Float64}
    z0::Vector{Float64}
    width::Vector{Float64}
    n::Vector{Float64}
    vacuum_density::Float64
    n_timesteps::Int
    tree::NearestNeighbors.KDTree
end

struct LineKDTree <: NNInterpolator
    r::Vector{Float64}
    z::Vector{Float64}
    n::Vector{Float64}
    width::Vector{Float64}
    line_tree::NearestNeighbors.KDTree
    line_tree_z::NearestNeighbors.KDTree
    n_timesteps::Int
end

r_min(grid::NNInterpolator) = minimum(grid.r)
r_max(grid::NNInterpolator) = maximum(grid.r)
z_min(grid::NNInterpolator) = minimum(grid.z)
z_max(grid::NNInterpolator) = maximum(grid.z)

function create_wind_kdtree(r, z, zmax, z0, width, density, vacuum_density, n_timesteps)
    points = hcat(r, z)'
    points = convert(Array{Float64,2}, points)
    kdtree = NearestNeighbors.KDTree(points)
    wind_kdtree =
        KDTree(r, z, zmax, z0, width, density, vacuum_density, n_timesteps, kdtree)
    return wind_kdtree
end


function create_wind_kdtree(
    integrators::Array{Sundials.IDAIntegrator},
    n_timesteps = 10000,
    vacuum_density = 1e2,
)
    r, z, zmax, z0, width, density =
        get_dense_solution_from_integrators(integrators, n_timesteps)
    return create_wind_kdtree(r, z, zmax, z0, width, density, vacuum_density, n_timesteps)
end

function create_line_kdtree(r, z, density, width, n_timesteps)
    points_line = hcat(r, z)'
    points_line = convert(Array{Float64,2}, points_line)
    line_tree = NearestNeighbors.KDTree(points_line)
    points_line = z'
    points_line = convert(Array{Float64,2}, points_line)
    line_tree_z = NearestNeighbors.KDTree(points_line)
    return LineKDTree(r, z, density, width, line_tree, line_tree_z, n_timesteps)
end

function create_line_kdtree(
    integrator::Sundials.IDAIntegrator,
    n_timesteps = 10000,
)
    r, z, zmax, z0, width, density =
        get_dense_solution_from_integrator(integrator, n_timesteps)
    return create_line_kdtree(r, z, density, width, n_timesteps)
end

function create_lines_kdtree(
    integrators::Array{Sundials.IDAIntegrator},
    n_timesteps = 10000,
    vacuum_density = 1e2,
)
    return [create_line_kdtree(integrator, n_timesteps) for integrator in integrators]
end

function get_closest_point(kdtree::KDTree, r, z)
    idcs, distance = NearestNeighbors.nn(kdtree.tree, [r, z])
    rc = kdtree.r[idcs]
    zc = kdtree.z[idcs]
    zmaxc = kdtree.zmax[idcs]
    z0c = kdtree.z0[idcs]
    widthc = kdtree.width[idcs]
    nc = kdtree.n[idcs]
    return rc, zc, zmaxc, z0c, widthc, nc, distance
end

function get_closest_points(kdtree::KDTree, points)
    idcs, distances = NearestNeighbors.nn(kdtree.tree, points)
    rc = kdtree.r[idcs]
    zc = kdtree.z[idcs]
    zmaxc = kdtree.zmax[idcs]
    z0c = kdtree.z0[idcs]
    widthc = kdtree.width[idcs]
    nc = kdtree.n[idcs]
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
Gets density of nearest calculated point. If the point is considered to be
outside the wind, then it returns the vacuum density.
"""
function get_density(kdtree::KDTree, r, z)
    rc, zc, zmax, z0, width, density, distance = get_closest_point(kdtree, r, z)
    if is_point_outside_wind(z, zmax, z0, distance, width)
        return kdtree.vacuum_density
    else
        return density
    end
end


function get_density_points(kdtree::KDTree, points)
    rc, zc, zmax, z0, width, density, distances = get_closest_points(kdtree, points)
    z = points[2, :]
    mask = is_point_outside_wind.(z, zmax, z0, distances, width) #(z .> zmax) || (distances .> width)
    density[mask] .= kdtree.vacuum_density
    density
end

function get_closest_density_distance(kdtree::KDTree, z)
    idcs, distance = NearestNeighbors.nn(kdtree.tree, [z])
    nc = kdtree.n[idcs]
    return nc, distance
end

#function get_density(
#    kdtree::KDTree,
#    kdtrees::Vector{KDTree},
#    kdtrees_z::Vector{KDTree},
#    r,
#    z,
#)
#    rc, zc, zmax, z0, width, density, distance = get_closest_point(kdtree, r, z)
#    if is_point_outside_wind(z, zmax, z0, distance, width)
#        return kdtree.vacuum_density
#    end
#    distances = []
#    densities = []
#    for (i, kdtree) in enumerate(kdtrees)
#        rc, zc, zmax, z0, width, density, distance = get_closest_point(kdtree, r, z)
#        if distance > width
#            continue
#        end
#        kdtree_z = kdtrees_z[i]
#        z_idx, _ = NearestNeighbors.nn(kdtree_z.tree, [z])
#        dens = kdtree_z.n[z_idx]
#        dist =
#            Distances.evaluate(Euclidean(), [r, z], [kdtree_z.r[z_idx], kdtree_z.z[z_idx]])
#        push!(densities, dens)
#        push!(distances, dist)
#    end
#    if length(distances) == 0
#        return 1e2
#    elseif length(distances) == 1
#        return densities[1]
#    else
#        idcs = sortperm(distances)[1:2]
#        d1, d2 = distances[idcs]
#        n1, n2 = densities[idcs]
#        return (n1 * d2 + n2 * d1) / (d1 + d2)
#    end
#end
#
