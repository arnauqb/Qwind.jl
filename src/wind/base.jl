using PyCall
import ConcaveHull, Interpolations, Sundials
export Wind, get_density

struct Wind{T,U,V}
    wind_hull::ConcaveHull.Hull
    density_grid::DensityGrid{T,U}
    velocity_grid::VelocityGrid{T,U}
    grid_iterator::GridIterator{T,U,V}
    vacuum_density::T
    update_grid_flag::UpdateGridFlag
    function Wind(wind_hull, density_grid, velocity_grid, vacuum_density, update_grid_flag)
        iterator = GridIterator(density_grid.r_range, density_grid.z_range)
        return new{Float64,Int,Bool}(
            wind_hull,
            density_grid,
            velocity_grid,
            iterator,
            vacuum_density,
            update_grid_flag,
        )
    end
end

Wind(nr, nz; vacuum_density = 1e2, update_grid_flag = AverageGrid()) = Wind(
    Hull(),
    DensityGrid(nr, nz, vacuum_density),
    VelocityGrid(nr, nz, 0.0),
    vacuum_density,
    update_grid_flag,
)

function Wind(parameters::Parameters)
    return Wind(
        parameters.radiation_grid_nr,
        parameters.radiation_grid_nz,
        vacuum_density = parameters.vacuum_density,
        update_grid_flag = parameters.update_grid_flag,
    )
end
function Wind(config::Dict)
    parameters = Parameters(config)
    return Wind(parameters)
end

function Wind(
    streamlines::Union{Streamlines,Nothing};
    nr = "auto",
    nz = 50,
    vacuum_density = 1e2,
    update_grid_flag = Average(),
)
    if nr != "auto"
        nr = Int(nr)
    end
    nz = Int(nz)
    if streamlines === nothing
        hull = Hull()
        density_grid = DensityGrid(nr, nz, vacuum_density)
        velocity_grid = VelocityGrid(nr, nz, 0.0)
    else
        hull = Hull(streamlines)
        velocity_grid = VelocityGrid(streamlines, hull, nr = nr, nz = nz)
        density_grid = DensityGrid(streamlines, hull, nr = nr, nz = nz)
    end
    return Wind(hull, density_grid, velocity_grid, vacuum_density, update_grid_flag)
end

function get_density(wi::Wind, r, z)
    if is_point_in_wind(wi, [r, z])
        return wi.grid.interpolator(r, z)
    else
        return wi.vacuum_density
    end
end
get_density(wi::Wind, point) = get_density(wi, point[1], point[2])

#function get_disk_density_from_streamlines(streamlines, r_range)
#    r0 = [line.r[1] for line in streamlines]
#    log_n0 = log10.([line.n[1] for line in streamlines])
#    println(log_n0)
#    interpolator = Interpolations.interpolate((r0,), log_n0, Gridded(Linear()))
#    interpolator = Interpolations.extrapolate(interpolator, 2)
#    n_interp = 10 .^ interpolator(r_range)
#    return n_interp
#end

function update_wind(model, streamlines::Streamlines)
    wind = model.wind
    @info "Updating wind grids... "
    flush()
    hull = Hull(streamlines)
    r, z, vr, vphi, vz = reduce_streamlines(streamlines)
    r0 = [line.r[1] for line in streamlines]
    r_range, z_range =
        get_spatial_grid(r, z, r0, wind.density_grid.nr, wind.density_grid.nz)
    disk_density = []
    for rp in r_range
        if rp < model.parameters.disk_r_out
            n0 = getn0(model, rp)
        else
            n0 = 1e2
        end
        push!(disk_density, n0)
    end
    velocity_grid =
        VelocityGrid(hull, r_range, z_range, r = r, z = z, vr = vr, vphi = vphi, vz = vz)
    density_grid = DensityGrid(
        hull,
        r_range,
        z_range,
        vr_grid = velocity_grid.vr_grid,
        vz_grid = velocity_grid.vz_grid,
        disk_density = disk_density,
    )
    new_wind =
        Wind(hull, density_grid, velocity_grid, wind.vacuum_density, wind.update_grid_flag)
    return new_wind
end

