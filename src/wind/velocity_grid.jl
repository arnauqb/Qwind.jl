using ConcaveHull
export VelocityGrid, get_velocity, interpolate_velocity, update_velocity_grid

struct VelocityGrid{T,U} <: InterpolationGrid
    r_range::Vector{T}
    z_range::Vector{T}
    vr_grid::Array{T,2}
    vphi_grid::Array{T,2}
    vz_grid::Array{T,2}
    nr::Union{U,String}
    nz::U
    vr_interpolator::Any
    vphi_interpolator::Any
    vz_interpolator::Any
    function VelocityGrid(
        r_range,
        z_range,
        vr_grid,
        vphi_grid,
        vz_grid,
        nr = nothing,
        nz = nothing,
    )
        if nr === nothing
            nr = length(r_range)
        end
        if nz === nothing
            nz = length(z_range)
        end
        vr_interpolator =
            Interpolations.interpolate((r_range, z_range), vr_grid, Gridded(Linear()))
        vphi_interpolator =
            Interpolations.interpolate((r_range, z_range), vphi_grid, Gridded(Linear()))
        vz_interpolator =
            Interpolations.interpolate((r_range, z_range), vz_grid, Gridded(Linear()))
        vr_interpolator = Interpolations.extrapolate(vr_interpolator, 0.0)
        vphi_interpolator = Interpolations.extrapolate(vphi_interpolator, 0.0)
        vz_interpolator = Interpolations.extrapolate(vz_interpolator, 0.0)
        return new{typeof(r_range[1]),Int}(
            r_range,
            z_range,
            vr_grid,
            vphi_grid,
            vz_grid,
            nr,
            nz,
            vr_interpolator,
            vphi_interpolator,
            vz_interpolator,
        )
    end
end

VelocityGrid(grid_data::Dict) = VelocityGrid(
    grid_data["r"],
    grid_data["z"],
    grid_data["vr_grid"],
    grid_data["vphi_grid"],
    grid_data["vz_grid"],
    get(grid_data, "nr", length(grid_data["r"])),
    #length(grid_data["r"]),
    get(grid_data, "nz", length(grid_data["z"])),
)

function VelocityGrid(h5_path, it_num)
    it_name = @sprintf "iteration_%03d" it_num
    grid_data = h5open(h5_path, "r") do file
        read(file, it_name * "/velocity_grid")
    end
    return VelocityGrid(grid_data)
end

function VelocityGrid(h5_path)
    it_keys = h5open(h5_path, "r") do file
        keys(read(file))
    end
    it_nums = [parse(Int, split(key, "_")[end]) for key in it_keys]
    return VelocityGrid(h5_path, maximum(it_nums))
end

function VelocityGrid(nr, nz, fill)
    r_range = [-1.0, 0.0]
    z_range = [-1.0, 0.0]
    vr_grid = zeros((2, 2))
    vphi_grid = zeros((2, 2))
    vz_grid = zeros((2, 2))
    return VelocityGrid(r_range, z_range, vr_grid, vphi_grid, vz_grid, nr, nz)
end

function VelocityGrid(
    hull::ConcaveHull.Hull,
    r,
    z,
    vr,
    vphi,
    vz,
    r0;
    nr = "auto",
    nz = 50,
    log = true,
)
    r_range, z_range = get_spatial_grid(r, z, r0, nr, nz, log = log)
    @info "Constructing velocity interpolators..."
    flush()
    vr_interp, vphi_interp, vz_interp = get_velocity_interpolators(r, z, vr, vphi, vz)
    @info "Done"
    @info "Filling velocity grids..."
    flush()
    r_grid = log10.(r_range)' .* ones(length(z_range))
    z_grid = log10.(z_range) .* ones(length(r_range))'
    vr_grid = vr_interp(r_grid, z_grid)
    vr_grid = reshape(vr_grid, length(z_range), length(r_range))'
    vz_grid = vz_interp(r_grid, z_grid)
    vz_grid = reshape(vz_grid, length(z_range), length(r_range))'
    vphi_grid = vphi_interp(r_grid, z_grid)
    vphi_grid = reshape(vphi_grid, length(z_range), length(r_range))'
    for (i, r) in enumerate(r_range)
        for (j, z) in enumerate(z_range)
            point = [r, z]
            if !is_point_in_wind(hull, point)
                vr_grid[i, j] = 0.0
                vphi_grid[i, j] = 0.0
                vz_grid[i, j] = 0.0
            end
        end
    end
    # add z = 0 line
    vr_grid = [vr_grid[:, 1] vr_grid]
    vphi_grid = [vphi_grid[:, 1] vphi_grid]
    vz_grid = [vz_grid[:, 1] vz_grid]
    pushfirst!(z_range, 0.0)
    grid = VelocityGrid(r_range, z_range, vr_grid, vphi_grid, vz_grid, nr, nz)
    @info "Done"
    return grid
end

function VelocityGrid(
    streamlines::Streamlines,
    hull::ConcaveHull.Hull;
    nr = "auto",
    nz = 50,
)
    r0 = [line.r[1] for line in streamlines]
    r, z, vr, vphi, vz, n = reduce_streamlines(streamlines)
    return VelocityGrid(
        hull,
        r,
        z,
        vr,
        vphi,
        vz,
        r0,
        nr = nr,
        nz = nz,
    )
end

function get_velocity(grid::VelocityGrid, r, z)
    if point_outside_grid(grid, r, z)
        return [0.0, 0.0, 0.0]
    end
    ridx = searchsorted_nearest(grid.r_range, r)
    zidx = searchsorted_nearest(grid.z_range, z)
    return [grid.vr_grid[ridx, zidx], grid.vphi_grid[ridx, zidx], grid.vz_grid[ridx, zidx]]
end

function interpolate_velocity(grid::VelocityGrid, r, z)
    if point_outside_grid(grid, r, z)
        return [0.0, 0.0, 0.0]
    end
    return [
        grid.vr_interpolator(r, z),
        grid.vphi_interpolator(r, z),
        grid.vz_interpolator(r, z),
    ]
end

function get_velocity_interpolators(r, z, vr, vphi, vz; type = "linear")
    mask = (r .> 0) .& (z .> 0)
    r = r[mask]
    z = z[mask]
    vr = vr[mask]
    vphi = vphi[mask]
    vz = vz[mask]
    r_log = log10.(r)
    z_log = log10.(z)
    points = hcat(r_log, z_log)
    vr_int = scipy_interpolate.LinearNDInterpolator(points, vr, fill_value = 0)
    vphi_int = scipy_interpolate.LinearNDInterpolator(points, vphi, fill_value = 0)
    vz_int = scipy_interpolate.LinearNDInterpolator(points, vz, fill_value = 0)
    return vr_int, vphi_int, vz_int
end


function update_velocity_grid(
    old_grid::VelocityGrid,
    update_method::UpdateGridFlag,
    streamlines::Streamlines,
    hull,
)
    return VelocityGrid(streamlines, hull, nr = old_grid.nr, nz = old_grid.nz)
end
