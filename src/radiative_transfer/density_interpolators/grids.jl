export get_spatial_grid,
    InterpolationGrid, get_density, get_density_points, get_density_grid

struct InterpolationGrid
    r_range::Vector{Float64}
    z_range::Vector{Float64}
    grid::Union{Array{Float64,2},Nothing}
    r_min::Float64
    z_min::Float64
    r_max::Float64
    z_max::Float64
    nr::Int
    nz::Int
    iterator::GridIterator
    function InterpolationGrid(r_range, z_range, grid)
        r_min = minimum(r_range)
        r_max = maximum(r_range)
        z_min = minimum(z_range)
        z_max = maximum(z_range)
        nr = length(r_range)
        nz = length(z_range)
        iterator = GridIterator(r_range, z_range)
        return new(r_range, z_range, grid, r_min, z_min, r_max, z_max, nr, nz, iterator)
    end
end

function is_point_outside_grid(grid::InterpolationGrid, r, z)
    if (z > grid.z_max) || (z < grid.z_min) || (r < grid.r_min) || (r > grid.r_max)
        return true
    else
        return false
    end
end

function get_spatial_grid(kdtree::KDTree, nr, nz)
    max_width = maximum(kdtree.width)
    r_min = kdtree.r[1] - kdtree.width[1] / 2
    #r_min = max(6, minimum(kdtree.r))
    r_max = maximum(kdtree.r) + max_width
    z_min = max(1e-6, minimum(kdtree.z))
    z_max = maximum(kdtree.z) + max_width
    r_range = 10 .^ range(log10(r_min), log10(r_max), length = nr)
    z_range = 10 .^ range(log10(z_min), log10(z_max), length = nz - 1)
    z_range = pushfirst!(z_range, 0.0)
    points = hcat([[r, z] for r in r_range for z in z_range]...)
    return r_range, z_range
end

function get_density(interpolator::GridInterpolator, r, z)
    return interpolator.interpolator(r, z)
    #if interpolator.grid.grid === nothing
    #    return interpolator.vacuum_density
    #end
    #if is_point_outside_grid(interpolator.grid, r, z)
    #    return interpolator.vacuum_density
    #end
    #r_idx = searchsorted_nearest(interpolator.grid.r_range, r)
    #z_idx = searchsorted_nearest(interpolator.grid.z_range, z)
    ##println("r_idx $r_idx $z_idx")
    #return interpolator.grid.grid[r_idx, z_idx]
end

function get_density_points(interpolator::GridInterpolator, points)
    return hcat([get_density(grid, point[1], point[2]) for point in eachcol(points)]...)
end
