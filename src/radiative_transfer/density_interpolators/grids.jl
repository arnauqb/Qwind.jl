import Interpolations
export get_spatial_grid,
    InterpolationGrid, get_density, get_density_points, get_density_grid

struct InterpolationGrid{T<:Float64}
    r_range::Vector{T}
    z_range::Vector{T}
    grid::Array{T,2}
    nr::Union{Int,String}
    nz::Int
    iterator::CellIterator{T}
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
        return new{typeof(r_range[1])}(
            r_range,
            z_range,
            grid,
            nr,
            nz,
            iterator,
            interpolator,
        )
    end
end

function point_outside_grid(grid::InterpolationGrid, r, z)
    if r < grid.r_range[1] || r > grid.r_range[end]
        return true
    end
    if z > grid.z_range[end]
        return true
    end
    return false
end

function get_density(grid::InterpolationGrid, r, z)
    if point_outside_grid(grid, r, z)
        return 1e2
    end
    ridx = searchsorted_nearest(grid.r_range, r)
    zidx = searchsorted_nearest(grid.z_range, z)
    return grid.grid[ridx, zidx]
end

get_density(grid::InterpolationGrid, point::Vector{Float64}) =
    get_density(grid, point[1], point[2])


