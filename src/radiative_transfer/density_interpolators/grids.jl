export get_spatial_grid,
    InterpolationGrid, get_density, get_density_points, get_density_grid

struct InterpolationGrid
    r_range::Union{Vector{Float64}, Nothing}
    z_range::Union{Vector{Float64}, Nothing}
    grid::Union{Array{Float64,2},Nothing}
    nr::Union{Int, String}
    nz::Int
    iterator::GridIterator
    function InterpolationGrid(r_range, z_range, grid, nr=nothing, nz=nothing)
        if nr === nothing
            nr = length(r_range)
        end
        if nz === nothing
            nz = length(z_range)
        end
        iterator = GridIterator(r_range, z_range)
        return new(r_range, z_range, grid, nr, nz, iterator)
    end
end

function get_spatial_grid(kdtree::KDTree, nr, nz)
    max_width = maximum(kdtree.width)
    r_min = kdtree.r[1] - kdtree.width[1] / 2
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
