#export VelocityGrid, get_velocity 
#
#struct VelocityGrid{T} <: InterpolationGrid{T}
#    r_range::Vector{T}
#    z_range::Vector{T}
#    grid_r::Array{T,2}
#    grid_z::Array{T,2}
#    nr::Union{Int,String}
#    nz::Int
#    iterator::CellIterator{T}
#    interpolator::Any
#    function VelocityGrid(r_range, z_range, grid_r, grid_z, nr = nothing, nz = nothing)
#        if nr === nothing
#            nr = length(r_range)
#        end
#        if nz === nothing
#            nz = length(z_range)
#        end
#        iterator = GridIterator(r_range, z_range)
#        if grid_r === nothing
#            interpolator = (r, z) -> 0.0
#        else
#            interpolator =
#                Interpolations.interpolate((r_range, z_range), grid, Gridded(Linear()))
#            interpolator = Interpolations.extrapolate(interpolator, 1e2)
#        end
#        return new{typeof(r_range[1])}(
#            r_range,
#            z_range,
#            grid,
#            nr,
#            nz,
#            iterator,
#            interpolator,
#        )
#    end
#end
#
#function get_density(grid::DensityGrid, r, z)
#    if point_outside_grid(grid, r, z)
#        return 1e2
#    end
#    ridx = searchsorted_nearest(grid.r_range, r)
#    zidx = searchsorted_nearest(grid.z_range, z)
#    return grid.grid[ridx, zidx]
#end
