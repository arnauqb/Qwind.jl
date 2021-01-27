using ImageFiltering
using Qwind

export SmoothedGrid, get_density, get_density_points, r_min, r_max, z_min, z_max

struct SmoothedGrid <: DensityInterpolator
    r_range::Vector{Float64}
    z_range::Vector{Float64}
    grid::Array{Float64,2}
    r_min::Float64
    z_min::Float64
    r_max::Float64
    z_max::Float64
    vacuum_density::Float64
    function SmoothedGrid(r_range, z_range, grid, vacuum_density = 1e2)
        r_min = minimum(r_range)
        r_max = maximum(r_range)
        z_min = minimum(z_range)
        z_max = maximum(z_range)
        grid = reshape(grid, (length(r_range), length(z_range)))
        return new(r_range, z_range, grid, r_min, z_min, r_max, z_max, vacuum_density)
    end
end

r_min(grid::SmoothedGrid) = grid.r_min
r_max(grid::SmoothedGrid) = grid.r_max
z_min(grid::SmoothedGrid) = grid.z_min
z_max(grid::SmoothedGrid) = grid.z_max

function SmoothedGrid(
    kdtree::Union{KDTree,Nothing};
    nr = 250,
    nz = 250,
    kernel_size = 10,
    vacuum_density = 1e2,
)
    if kdtree === nothing
        return SmoothedGrid([0.0], [0.0], [0.0], vacuum_density)
    end
    max_width = maximum(kdtree.width)
    r_min = max(6, minimum(kdtree.r) - max_width)
    r_max = maximum(kdtree.r) + max_width
    z_min = max(1e-6, minimum(kdtree.z) - max_width)
    z_max = maximum(kdtree.z) + max_width
    r_range = 10 .^ range(log10(r_min), log10(r_max), length = nr)
    z_range = 10 .^ range(log10(z_min), log10(z_max), length = nz)
    points = hcat([[r, z] for r in r_range for z in z_range]...)
    density_grid = reshape(get_density_points(kdtree, points), (nz, nr))'
    density_grid = imfilter(density_grid, Kernel.gaussian(kernel_size))
    return SmoothedGrid(r_range, z_range, density_grid, vacuum_density)
end

function is_point_outside_grid(grid::SmoothedGrid, r, z)
    if (z > grid.z_max) || (z < grid.z_min) || (r < grid.r_min) || (r > grid.r_max)
        return true
    else
        return false
    end
end

function get_density(grid::SmoothedGrid, r, z)
    if is_point_outside_grid(grid, r, z)
        return grid.vacuum_density
    end
    r_idx = searchsorted_nearest(grid.r_range, r)
    z_idx = searchsorted_nearest(grid.z_range, z)
    return grid.grid[r_idx, z_idx]
end

function get_density_points(grid::SmoothedGrid, points)
    return hcat([get_density(grid, point[1], point[2]) for point in eachcol(points)]...)
end

