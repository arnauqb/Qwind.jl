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
        return new{typeof(r_range[1]),Int}(
            r_range,
            z_range,
            grid,
            nr,
            nz,
            interpolator,
        )
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
    hull::ConcaveHull.Hull,
    r_range,
    z_range;
    vr_grid,
    vz_grid,
    disk_density,
)
    @info "Filling density grid..."
    flush()
    r_grid = log10.(r_range)' .* ones(length(z_range))
    z_grid = log10.(z_range) .* ones(length(r_range))'
    density_grid = 1e2 .* ones(length(r_range), length(z_range))
    density_grid[:,1] = disk_density
    for i = 2:length(r_range)
        ri = r_range[i]
        rii = r_range[i - 1]
        for j = 2:length(z_range)
            zi = z_range[j]
            zii = z_range[j - 1]
            point = [ri, zi]
            if !is_point_in_wind(hull, point)
                density_grid[i, j] = 1e2
            else
                D =
                    (ri * vr_grid[i, j]) / (2 * (ri^2 - rii^2)) +
                    (vz_grid[i, j]) / (zi - zii)
                rho =
                    (rii * vr_grid[i - 1, j] * density_grid[i - 1, j]) /
                    (2 * (ri^2 - rii^2))
                rho = rho + (vz_grid[i, j - 1] * density_grid[i, j - 1]) / (zi - zii)
                rho = abs(rho / D)
                if isnan(rho) || isinf(rho)
                    rho = 1e2
                end
                #rho = max(1e2, rho)
                density_grid[i, j] = rho
            end
        end
    end
    # add z = 0 line
    density_grid = [density_grid[:, 1] density_grid]
    z_save = copy(z_range)
    pushfirst!(z_save, 0.0)
    nr = length(r_range)
    nz = length(z_range)
    grid = DensityGrid(r_range, z_save, density_grid, nr, nz)
    @info "Done"
    flush()
    return grid
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
