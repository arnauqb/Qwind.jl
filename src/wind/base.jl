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
        density_grid = DensityGrid(streamlines, hull, nr = nr, nz = nz)
        velocity_grid = VelocityGrid(streamlines, hull, nr = nr, nz = nz)
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

function update_wind(wind::Wind, streamlines::Streamlines, integrators)
    @info "Updating wind grids... "
    flush()
    if maximum(wind.density_grid.grid) == wind.vacuum_density
        # first iteration, do not average
        new_wind = Wind(
            streamlines,
            nr = wind.density_grid.nr,
            nz = wind.density_grid.nz,
            vacuum_density = wind.vacuum_density,
            update_grid_flag = wind.update_grid_flag,
        )
        return new_wind
    end
    hull = Hull(streamlines)
    density_grid = update_density_grid(
        wind.density_grid,
        wind.update_grid_flag,
        streamlines,
        integrators,
        hull,
    )
    velocity_grid =
        update_velocity_grid(wind.velocity_grid, wind.update_grid_flag, streamlines, hull)
    new_wind =
        Wind(hull, density_grid, velocity_grid, wind.vacuum_density, wind.update_grid_flag)
    return new_wind
end

