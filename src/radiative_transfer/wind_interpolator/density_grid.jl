import ConcaveHull
export DensityGrid, get_density, update_density_grid, interpolate_density

struct DensityGrid{T} <: InterpolationGrid{T}
    r_range::Vector{T}
    z_range::Vector{T}
    grid::Array{T,2}
    nr::Union{Int,String}
    nz::Int
    iterator::CellIterator{T}
    interpolator::Any
    function DensityGrid(r_range, z_range, grid, nr = nothing, nz = nothing)
        if nr === nothing
            nr = length(r_range)
        end
        if nz === nothing
            nz = length(z_range)
        end
        iterator = GridIterator(r_range, z_range)
        interpolator =
            Interpolations.interpolate((r_range, z_range), grid, Gridded(Linear()))
        interpolator = Interpolations.extrapolate(interpolator, 1e2)
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

function get_density(grid::DensityGrid, r, z)
    if point_outside_grid(grid, r, z)
        return 1e2
    end
    ridx = searchsorted_nearest(grid.r_range, r)
    zidx = searchsorted_nearest(grid.z_range, z)
    return grid.grid[ridx, zidx]
end

function interpolate_density(grid::DensityGrid, r, z)
    if point_outside_grid(grid, r, z)
        return 1e2
    end
    return grid.interpolator(r,z)
end

get_density(grid::DensityGrid, point::Vector{Float64}) =
    get_density(grid, point[1], point[2])

function get_density_interpolator(
    r::Vector{Float64},
    z::Vector{Float64},
    n::Vector{Float64};
    type = "linear"
)
    mask = (r .> 0) .& (z .> 0)
    r = r[mask]
    z = z[mask]
    n = n[mask]
    r_log = log10.(r)
    z_log = log10.(z)
    log_n = log10.(n)
    points = hcat(r_log, z_log)
    if type == "linear"
        interp = scipy_interpolate.LinearNDInterpolator(points, log_n, fill_value = 2)
    elseif type == "nn"
        interp = scipy_interpolate.NearestNDInterpolator(points, log_n, rescale=true)
    elseif type == "clough_tocher"
        interp = scipy_interpolate.CloughTocher2DInterpolator(points, log_n, fill_value=2)
    elseif type == "rbf"
        interp = scipy_interpolate.Rbf(r_log, z_log, log_n)
    else
        error("interpolation type $type not supported")
    end
    return interp 
end

function construct_density_grid(nr, nz, vacuum_density = 1e2)
    r_range = zeros(2)
    z_range = zeros(2)
    density_grid = vacuum_density .* [[1.0, 1.0] [1.0, 1.0]]
    return DensityGrid(r_range, z_range, density_grid, nr, nz)
end

function construct_density_grid(
    r::Vector{Float64},
    z::Vector{Float64},
    n::Vector{Float64},
    r0s::Vector{Float64},
    hull::ConcaveHull.Hull;
    nr = "auto",
    nz = 50,
    log=true,
    interpolation_type = "linear",
)
    @info "Constructing density interpolator..."
    flush()
    interpolator = get_density_interpolator(r, z, n, type=interpolation_type)
    @info "Done"
    @info "Filling density grid..."
    flush()
    r_range, z_range = get_spatial_grid(r, z, r0s, nr, nz, log=log)
    r_range_log = log10.(r_range)
    z_range_log = log10.(z_range)
    density_grid = 1e2 .* ones((length(r_range), length(z_range)))
    for (i, r) in enumerate(r_range)
        for (j, z) in enumerate(z_range)
            point = [r, z]
            if !is_point_in_wind(hull, point)
                density_grid[i, j] = 1e2
            else
                density_grid[i, j] = 10 .^ interpolator(log10(r), log10(z))[1]
            end
        end
    end
    # add z = 0 line
    density_grid = [density_grid[:, 1] density_grid]
    pushfirst!(z_range, 0.0)
    grid = DensityGrid(r_range, z_range, density_grid, nr, nz)
    @info "Done"
    flush()
    return grid
end

function update_density_grid(
    old_grid::DensityGrid,
    hull::ConcaveHull.Hull,
    r::Vector{Float64},
    z::Vector{Float64},
    n::Vector{Float64},
    r0s::Vector{Float64},
)
    r_range, z_range = get_spatial_grid(r, z, r0s, old_grid.nr, old_grid.nz)
    r_range_log = log10.(r_range)
    z_range_log = log10.(z_range)
    interpolator = get_density_interpolator(r, z, n)
    density_grid = 1e2 .* ones((length(r_range), length(z_range)))
    @info "Averaging grids..."
    flush()
    for (i, r) in enumerate(r_range)
        for (j, z) in enumerate(z_range)
            point = [r, z]
            if !is_point_in_wind(hull, point)
                density_grid[i, j] = 1e2
            else
                density_grid[i, j] =
                    10 .^ (
                        (
                            interpolator(log10(r), log10(z))[1] +
                            log10(get_density(old_grid, r, z))
                        ) / 2.0
                    )
            end
        end
    end
    # add z = 0 line
    density_grid = [density_grid[:, 1] density_grid]
    pushfirst!(z_range, 0.0)
    grid = DensityGrid(r_range, z_range, density_grid, old_grid.nr, old_grid.nz)
    return grid
end
