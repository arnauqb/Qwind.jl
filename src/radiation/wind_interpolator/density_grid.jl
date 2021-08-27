using HDF5, Printf
import ConcaveHull
export DensityGrid, get_density, update_density_grid, interpolate_density

struct DensityGrid{T, U, V} <: InterpolationGrid
    r_range::Vector{T}
    z_range::Vector{T}
    grid::Array{T,2}
    nr::Union{U, String}
    nz::U
    iterator::GridIterator{T, U, V}
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
        return new{typeof(r_range[1]), Int, Bool}(
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

DensityGrid(grid_data::Dict) =
    DensityGrid(grid_data["r"], grid_data["z"], grid_data["grid"])

function DensityGrid(h5_path::String, it_num)
    it_name = @sprintf "iteration_%03d" it_num
    grid_data = h5open(h5_path, "r") do file
        read(file, it_name * "/density_grid")
    end
    return DensityGrid(grid_data)
end

function DensityGrid(h5_path::String)
    it_keys = h5open(h5_path, "r") do file
        keys(read(file))
    end
    it_nums = [parse(Int, split(key, "_")[end]) for key in it_keys]
    return DensityGrid(h5_path, maximum(it_nums))
end


function DensityGrid(nr::Union{String,Int}, nz::Int, vacuum_density::Float64)
    r_range = [-1.0, 0.0]
    z_range = [-1.0, 0.0]
    density_grid = vacuum_density .* [[1.0, 1.0] [1.0, 1.0]]
    return DensityGrid(r_range, z_range, density_grid, nr, nz)
end

function DensityGrid(
    r::Vector{Float64},
    z::Vector{Float64},
    n::Vector{Float64},
    r0s::Vector{Float64},
    hull::ConcaveHull.Hull;
    nr = "auto",
    nz = 50,
    log = true,
    interpolation_type = "linear",
)
    @info "Constructing density interpolator..."
    flush()
    interpolator = get_density_interpolator(r, z, n, type = interpolation_type)
    @info "Done"
    @info "Filling density grid..."
    flush()
    r_range, z_range = get_spatial_grid(r, z, r0s, nr, nz, log = log)
    r_grid = log10.(r_range)' .* ones(length(z_range))
    z_grid = log10.(z_range) .* ones(length(r_range))'
    density_grid = interpolator(r_grid, z_grid)
    density_grid = 10 .^ reshape(density_grid, length(z_range), length(r_range))'
    for (i, r) in enumerate(r_range)
        for (j, z) in enumerate(z_range)
            point = [r, z]
            if !is_point_in_wind(hull, point)
                density_grid[i, j] = 1e2
            end
        end
    end
    # add z = 0 line
    #density_grid = [density_grid[:, 1] density_grid]
    #pushfirst!(z_range, 0.0)
    grid = DensityGrid(r_range, z_range, density_grid, nr, nz)
    @info "Done"
    flush()
    return grid
end

function DensityGrid(
    integrators::Vector{<:Sundials.IDAIntegrator},
    max_times,
    hull;
    nr = "auto",
    nz = 50,
    interpolation_type = "linear",
)
    r0 = [integ.p.r0 for integ in integrators]
    integrators_interpolated = interpolate_integrators(
        integrators,
        max_times = max_times,
        n_timesteps = 1000,
        log = true,
    )
    r, z, vr, vphi, vz, n = reduce_integrators(integrators_interpolated)
    return DensityGrid(
        r,
        z,
        n,
        r0,
        hull,
        nr = nr,
        nz = nz,
        log = true,
        interpolation_type = "linear",
    )
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

function interpolate_density(grid::DensityGrid, r, z)
    if point_outside_grid(grid, r, z)
        return 1e2
    end
    return grid.interpolator(r, z)
end


function get_density_interpolator(
    r::Vector{Float64},
    z::Vector{Float64},
    n::Vector{Float64};
    type = "linear",
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
        interp = scipy_interpolate.NearestNDInterpolator(points, log_n, rescale = true)
    elseif type == "clough_tocher"
        interp = scipy_interpolate.CloughTocher2DInterpolator(points, log_n, fill_value = 2)
    elseif type == "rbf"
        interp = scipy_interpolate.Rbf(r_log, z_log, log_n)
    else
        error("interpolation type $type not supported")
    end
    return interp
end

function update_density_grid(
    old_grid::DensityGrid,
    update_method::AverageGrid,
    hull::ConcaveHull.Hull,
    r::Vector{Float64},
    z::Vector{Float64},
    n::Vector{Float64},
    r0::Vector{Float64},
)
    r_range, z_range = get_spatial_grid(r, z, r0, old_grid.nr, old_grid.nz)
    interpolator = get_density_interpolator(r, z, n)
    density_grid = 1e2 .* ones((length(r_range), length(z_range)))
    @info "Averaging grids..."
    flush()
    r_grid = log10.(r_range)' .* ones(length(z_range))
    z_grid = log10.(z_range) .* ones(length(r_range))'
    density_grid = interpolator(r_grid, z_grid)
    density_grid = 10 .^ reshape(density_grid, length(z_range), length(r_range))'
    for (i, r) in enumerate(r_range)
        for (j, z) in enumerate(z_range)
            point = [r, z]
            if !is_point_in_wind(hull, point)
                density_grid[i, j] = 
                    10 .^ (
                        (
                            log10(1e2) +
                            log10(get_density(old_grid, r, z))
                        ) / 2.0
                    )
            else
                density_grid[i, j] = 
                    10 .^ (
                        (
                            log10(density_grid[i,j]) +
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

function update_density_grid(
    old_grid::DensityGrid,
    update_method::ReplaceGrid,
    hull::ConcaveHull.Hull,
    r::Vector{Float64},
    z::Vector{Float64},
    n::Vector{Float64},
    r0::Vector{Float64},
)
    r_range, z_range = get_spatial_grid(r, z, r0, old_grid.nr, old_grid.nz)
    interpolator = get_density_interpolator(r, z, n)
    density_grid = 1e2 .* ones((length(r_range), length(z_range)))
    @info "Replacing grids..."
    flush()
    r_grid = log10.(r_range)' .* ones(length(z_range))
    z_grid = log10.(z_range) .* ones(length(r_range))'
    density_grid = interpolator(r_grid, z_grid)
    density_grid = 10 .^ reshape(density_grid, length(z_range), length(r_range))'
    for (i, r) in enumerate(r_range)
        for (j, z) in enumerate(z_range)
            point = [r, z]
            if !is_point_in_wind(hull, point)
                density_grid[i, j] = 1e2
            end
        end
    end
    # add z = 0 line
    density_grid = [density_grid[:, 1] density_grid]
    pushfirst!(z_range, 0.0)
    grid = DensityGrid(r_range, z_range, density_grid, old_grid.nr, old_grid.nz)
    return grid
end

function update_density_grid(
    old_grid::DensityGrid,
    update_method::UpdateGridFlag,
    integrators::Vector{<:Sundials.IDAIntegrator},
    max_times,
    hull,
)
    r0 = [integ.p.r0 for integ in integrators]
    integrators_interpolated = interpolate_integrators(
        integrators,
        max_times = max_times,
        n_timesteps = 1000,
        log = true,
    )
    r, z, vr, vphi, vz, n = reduce_integrators(integrators_interpolated)
    return update_density_grid(old_grid, update_method, hull, r, z, n, r0)
end
