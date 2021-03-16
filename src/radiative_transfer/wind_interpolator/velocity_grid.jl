export VelocityGrid, get_velocity, interpolate_velocity

struct VelocityGrid{T} <: InterpolationGrid{T}
    r_range::Vector{T}
    z_range::Vector{T}
    vr_grid::Array{T,2}
    vz_grid::Array{T,2}
    nr::Union{Int,String}
    nz::Int
    iterator::CellIterator{T}
    vr_interpolator::Any
    vz_interpolator::Any
    function VelocityGrid(r_range, z_range, vr_grid, vz_grid, nr = nothing, nz = nothing)
        if nr === nothing
            nr = length(r_range)
        end
        if nz === nothing
            nz = length(z_range)
        end
        iterator = GridIterator(r_range, z_range)
        vr_interpolator =
            Interpolations.interpolate((r_range, z_range), vr_grid, Gridded(Linear()))
        vz_interpolator =
            Interpolations.interpolate((r_range, z_range), vz_grid, Gridded(Linear()))
        vr_interpolator = Interpolations.extrapolate(vr_interpolator, 0.0)
        vz_interpolator = Interpolations.extrapolate(vz_interpolator, 0.0)
        return new{typeof(r_range[1])}(
            r_range,
            z_range,
            vr_grid,
            vz_grid,
            nr,
            nz,
            iterator,
            vr_interpolator,
            vz_interpolator
        )
    end
end

function get_velocity(grid::VelocityGrid, r, z)
    if point_outside_grid(grid, r, z)
        return [0.0, 0.0]
    end
    ridx = searchsorted_nearest(grid.r_range, r)
    zidx = searchsorted_nearest(grid.z_range, z)
    return [grid.vr_grid[ridx, zidx], grid.vz_grid[ridx, zidx]]
end

function interpolate_velocity(grid::VelocityGrid, r, z)
    if point_outside_grid(grid, r, z)
        return [0.0, 0.0]
    end
    return [grid.vr_interpolator(r,z), grid.vz_interpolator(r,z)]
end

function get_velocity_interpolators(
    r::Vector{Float64},
    z::Vector{Float64},
    vr::Vector{Float64},
    vz::Vector{Float64};
    type = "linear"
)
    r_log = log10.(r)
    z_log = log10.(z)
    points = hcat(r_log, z_log)
    if type == "linear"
        vr_int = scipy_interpolate.LinearNDInterpolator(points, vr, fill_value = 0)
        vz_int = scipy_interpolate.LinearNDInterpolator(points, vz, fill_value = 0)
    elseif type == "nn"
        vr_int = scipy_interpolate.NearestNDInterpolator(points, vr)
        vz_int = scipy_interpolate.NearestNDInterpolator(points, vz)
    else
        error("interpolation type $type not supported")
    end
    return vr_int, vz_int
end

function construct_velocity_grid(nr, nz)
    r_range = zeros(2)
    z_range = zeros(2)
    vr_grid = zeros((2,2))
    vz_grid = zeros((2,2))
    return VelocityGrid(r_range, z_range, vr_grid, vz_grid, nr, nz)
end

function construct_velocity_grid(
    r::Vector{Float64},
    z::Vector{Float64},
    vr::Vector{Float64},
    vz::Vector{Float64},
    r0s::Vector{Float64},
    hull::ConcaveHull.Hull;
    nr = "auto",
    nz = 50,
    log=true,
    interpolation_type="linear"
)
    r_range, z_range = get_spatial_grid(r, z, r0s, nr, nz, log=log)
    @info "Constructing velocity interpolators..."
    flush()
    vr_interp, vz_interp = get_velocity_interpolators(r, z, vr, vz, type=interpolation_type)
    @info "Done"
    @info "Filling velocity grids..."
    flush()
    r_range_log = log10.(r_range)
    z_range_log = log10.(z_range)
    vr_grid = zeros((length(r_range), length(z_range)))
    vz_grid = zeros((length(r_range), length(z_range)))
    for (i, r) in enumerate(r_range)
        for (j, z) in enumerate(z_range)
            point = [r, z]
            if !is_point_in_wind(hull, point)
                continue
            else
                vr_grid[i,j] = vr_interp(log10(r), log10(z))[1]
                vz_grid[i,j] = vz_interp(log10(r), log10(z))[1]
            end
        end
    end
    # add z = 0 line
    vr_grid = [vr_grid[:,1] vr_grid]
    vz_grid = [vz_grid[:,1] vz_grid]
    pushfirst!(z_range, 0.0)
    grid = VelocityGrid(r_range, z_range, vr_grid, vz_grid, nr, nz)
    @info "Done"
    return grid
end

