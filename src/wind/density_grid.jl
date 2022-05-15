using HDF5, Printf
import ConcaveHull
export DensityGrid, get_density, update_density_grid, interpolate_density

struct DensityGrid{T,U} <: InterpolationGrid
    r_range::Vector{T}
    z_range::Vector{T}
    grid::Array{T,2}
    nr::Union{U,String}
    nz::U
    interpolator::Any
    function DensityGrid(r_range, z_range, grid, nr = nothing, nz = nothing)
        if nr === nothing
            nr = length(r_range)
        end
        if nz === nothing
            nz = length(z_range)
        end
        interpolator =
            Interpolations.interpolate((r_range, z_range), grid, Gridded(Linear()))
        interpolator = Interpolations.extrapolate(interpolator, 1e2)
        return new{typeof(r_range[1]),Int}(r_range, z_range, grid, nr, nz, interpolator)
    end
end

DensityGrid(grid_data::Dict) = DensityGrid(
    grid_data["r"],
    grid_data["z"],
    grid_data["grid"],
    get(grid_data, "nr", length(grid_data["r"])),
    get(grid_data, "nz", length(grid_data["z"])),
)

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
)
    @info "Constructing density interpolator..."
    flush()
    interpolator = get_density_interpolator(r, z, n)
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
    density_grid = [density_grid[:, 1] density_grid]
    pushfirst!(z_range, 0.0)
    grid = DensityGrid(r_range, z_range, density_grid, nr, nz)
    @info "Done"
    flush()
    return grid
end

function DensityGrid(streamlines, hull; nr = "auto", nz = 50)
    r0 = [line.r[1] for line in streamlines]
    r, z, vr, vphi, vz, n = reduce_streamlines(streamlines)
    return DensityGrid(r, z, n, r0, hull, nr = nr, nz = nz, log = true)
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


function get_density_interpolator(r, z, n)
    mask = (r .> 0) .& (z .> 0)
    n = n[mask]
    r = r[mask]
    z = z[mask]
    r_log = log10.(r)
    z_log = log10.(z)
    log_n = log10.(n)
    points = hcat(r_log, z_log)
    interp = scipy_interpolate.LinearNDInterpolator(points, log_n, fill_value = 2)
    return interp
end

function update_density_grid(
    old_grid::DensityGrid,
    update_method::AverageGrid,
    hull::ConcaveHull.Hull,
    r::Vector{Float64},
    z::Vector{Float64},
    r0::Vector{Float64},
    n,
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
                    10 .^ ((log10(1e2) + log10(get_density(old_grid, r, z))) / 2.0)
            else
                density_grid[i, j] =
                    10 .^
                    ((log10(density_grid[i, j]) + log10(get_density(old_grid, r, z))) / 2.0)
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
    r0::Vector{Float64},
    n,
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

function calculate_densities(streamlines)
    A = zeros(2, 2)
    b = zeros(2)
    # get line widths
    lw_pairs = get_width_pairs(streamlines[2], streamlines[1], streamlines[3], A, b)
    lws = [[[Inf, lw_pairs[i][1]] for i = 1:length(lw_pairs)]]
    for i = 2:(length(streamlines) - 1)
        lw_pairs =
            get_width_pairs(streamlines[i], streamlines[i - 1], streamlines[i + 1], A, b)
        push!(lws, lw_pairs)
    end
    lws_last = [[lws[end][i][2], Inf] for i = 1:length(lws[end])]
    push!(lws, lws_last)
    # calculate densities
    rs = Float64[]
    zs = Float64[]
    ns = Float64[]
    for (i, sl) in enumerate(streamlines)
        r0 = sl.r[1]
        n0 = sl.n[1]
        v0 = sl.vz[1]
        lw0 = sl.width[1]
        for j = 1:(length(sl.r) - 1)
            r = sl.r[j]
            z = sl.z[j]
            d = sqrt(r^2 + z^2)
            v = sqrt(sl.vr[j]^2 + sl.vz[j]^2)
            left_lw = min(lws[i][j][1], lw0)
            right_lw = min(lws[i][j][2], lw0)
            lw = left_lw + right_lw
            n = (r0 * lw0 * v0 * n0) / (d * v * lw)
            push!(ns, n)
            push!(rs, r)
            push!(zs, z)
        end
    end
    return rs, zs, ns
end

function update_density_grid(
    old_grid::DensityGrid,
    update_method::UpdateGridFlag,
    streamlines::Streamlines,
    hull,
)
    #r, z, _, _, _, _ = reduce_streamlines(streamlines)
    r, z, n = calculate_densities(streamlines)
    r0s = [line.r[1] for line in streamlines]
    return update_density_grid(
        old_grid,
        update_method,
        hull,
        r,
        z,
        r0s,
        n,
    )
end
