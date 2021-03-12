using PyCall
import ConcaveHull
export DensityGrid, get_density, update_density_grid

const scipy_interpolate = PyNULL()

function __init__()
    copy!(scipy_interpolate, pyimport_conda("scipy.interpolate", "scipy"))
end

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

function get_density(grid::DensityGrid, r, z)
    if point_outside_grid(grid, r, z)
        return 1e2
    end
    ridx = searchsorted_nearest(grid.r_range, r)
    zidx = searchsorted_nearest(grid.z_range, z)
    return grid.grid[ridx, zidx]
end

get_density(grid::DensityGrid, point::Vector{Float64}) =
    get_density(grid, point[1], point[2])

function get_density_interpolator(
    r::Vector{Float64},
    z::Vector{Float64},
    n::Vector{Float64},
)
    r_log = log10.(r)
    z_log = log10.(z)
    log_n = log10.(n)
    points = hcat(r_log, z_log)
    linear_int = scipy_interpolate.LinearNDInterpolator(points, log_n, fill_value = 2)
    return linear_int
end

function construct_density_grid(nr, nz, vacuum_density=1e2)
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
)
    interpolator = get_density_interpolator(r, z, n)
    r_range, z_range = get_spatial_grid(r, z, r0s, nr, nz)
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
    return grid
end

function get_spatial_grid(
    r::Vector{Float64},
    z::Vector{Float64},
    r0s::Vector{Float64},
    nr = "auto",
    nz = 50,
)
    r_min = max(6.0, minimum(r))
    z_min = max(minimum(z), 1e-6)
    r_max = min(1e4, maximum(r))
    z_max = min(1e4, maximum(z))
    if nr == "auto"
        r_range = r0s
        additional_r = collect(range(r_range[end], r_max, step = 100)[2:end])
        r_range = vcat(r_range, additional_r)
    else
        r_range = 10 .^ range(log10(r_min), log10(r_max), length = nr)
    end
    z_range = 10 .^ range(log10(z_min), log10(z_max), length = nz - 1)
    r_range = round.(r_range, digits = 7)
    z_range = round.(z_range, digits = 7)
    points = hcat([[r, z] for r in r_range for z in z_range]...)
    return r_range, z_range
end

function update_density_grid(old_grid::DensityGrid, hull::ConcaveHull.Hull, integrators)
    r0s = [integ.p.r0 for integ in integrators]
    r, z, n = reduce_integrators(integrators, n_timesteps = 1000)
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
