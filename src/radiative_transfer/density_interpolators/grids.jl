import Interpolations
export get_spatial_grid,
    InterpolationGrid, get_density, get_density_points, get_density_grid

struct InterpolationGrid
    r_range::Union{Vector{Float64},Nothing}
    z_range::Union{Vector{Float64},Nothing}
    grid::Union{Array{Float64,2},Nothing}
    nr::Union{Int,String}
    nz::Int
    iterator::CellIterator
    interpolator::Any
    function InterpolationGrid(r_range, z_range, grid, nr = nothing, nz = nothing)
        if nr === nothing
            nr = length(r_range)
        end
        if nz === nothing
            nz = length(z_range)
        end
        iterator = GridIterator(r_range, z_range)
        if grid === nothing
            interpolator = (r, z) -> 1e2
        else
            interpolator =
                Interpolations.interpolate((r_range, z_range), grid, Gridded(Linear()))
            interpolator = Interpolations.extrapolate(interpolator, 1e2)
        end
        return new(r_range, z_range, grid, nr, nz, iterator, interpolator)
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

function get_density(grid::InterpolationGrid, r, z)
    return grid.interpolator(r, z)
end

get_density(grid::InterpolationGrid, point::Vector{Float64}) =
    get_density(grid, point[1], point[2])


