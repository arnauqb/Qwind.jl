import NearestNeighbors
using Distances, Sundials
export create_wind_kdtree,
    get_density,
    r_min,
    r_max,
    z_min,
    z_max,
    KDTree,
    get_closest_points,
    get_closest_point,
    get_density_points,
    is_point_outside_wind

struct KDTree <: DensityInterpolator
    r::Vector{Float64}
    z::Vector{Float64}
    zmax::Vector{Float64}
    z0::Vector{Float64}
    width::Vector{Float64}
    n::Vector{Float64}
    vacuum_density::Float64
    tree::NearestNeighbors.KDTree
end

r_min(grid::KDTree) = minimum(grid.r) 
r_max(grid::KDTree) = maximum(grid.r) 
z_min(grid::KDTree) = minimum(grid.z) 
z_max(grid::KDTree) = maximum(grid.z) 

function create_wind_kdtree(r, z, zmax, z0, width, density, vacuum_density = 1e2)
    points = hcat(r, z)'
    points = convert(Array{Float64,2}, points)
    kdtree = NearestNeighbors.KDTree(points)
    wind_kdtree = KDTree(r, z, zmax, z0, width, density, vacuum_density, kdtree)
    return wind_kdtree
end

function create_wind_kdtree(
    integrators::Array{Sundials.IDAIntegrator},
    n_timesteps = 10000,
    vacuum_density = 1e2,
)
    r, z, zmax, z0, width, density =
        get_dense_solution_from_integrators(integrators, n_timesteps)
    return create_wind_kdtree(r, z, zmax, z0, width, density, vacuum_density)
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
