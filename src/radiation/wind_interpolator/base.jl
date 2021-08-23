using PyCall
import ConcaveHull, Interpolations, Sundials
export WindInterpolator, get_density

struct WindInterpolator{T, U, V}
    wind_hull::ConcaveHull.Hull
    density_grid::DensityGrid{T, U, V}
    velocity_grid::VelocityGrid{T, U, V}
    vacuum_density::T
    n_timesteps::Int
    update_grid_flag::UpdateGridFlag
end

WindInterpolator(nr, nz; vacuum_density = 1e2, n_timesteps = 1000, update_grid_flag = AverageGrid()) = WindInterpolator(
    Hull(),
    DensityGrid(nr, nz, vacuum_density),
    VelocityGrid(nr, nz, 0.0),
    vacuum_density,
    n_timesteps,
    update_grid_flag
)
function WindInterpolator(config::Dict)
    update_method = get(config, :update_grid_method, "average")
    if update_method == "average"
        update_method = AverageGrid()
    elseif update_method == "replace"
        update_method = ReplaceGrid()
    else
        error("Grid update method $update_method not recoginsed")
    end
    return WindInterpolator(
        config[:nr],
        config[:nz],
        vacuum_density = config[:vacuum_density],
        n_timesteps = config[:n_integrator_interpolation],
        update_grid_flag = update_method
    )
end

function WindInterpolator(
    integrators;
    nr = "auto",
    nz = 50,
    vacuum_density = 1e2,
    n_timesteps = 1000,
    update_grid_flag = Average(),
)
    if nr != "auto"
        nr = Int(nr)
    end
    nz = Int(nz)
    if integrators === nothing
        hull = Hull()
        density_grid = DensityGrid(nr, nz, vacuum_density)
        velocity_grid = VelocityGrid(nr, nz, 0.0)
    else
        r0 = [integ.p.r0 for integ in integrators]
        max_times = get_intersection_times(integrators)
        hull = Hull(integrators, max_times)
        density_grid = DensityGrid(integrators, max_times, hull, nr = nr, nz = nz)
        velocity_grid = VelocityGrid(integrators, max_times, hull, nr = nr, nz = nz)
    end
    return WindInterpolator(hull, density_grid, velocity_grid, vacuum_density, n_timesteps, update_grid_flag)
end


function update_wind_interpolator(wi::WindInterpolator, integrators::Vector{<:Sundials.IDAIntegrator})
    if maximum(wi.density_grid.grid) == wi.vacuum_density
        # first iteration, do not average
        return WindInterpolator(
            integrators,
            nr = wi.density_grid.nr,
            nz = wi.density_grid.nz,
            vacuum_density = wi.vacuum_density,
            n_timesteps = wi.n_timesteps,
            update_grid_flag = wi.update_grid_flag,
        )
    end
    r0 = [integ.p.r0 for integ in integrators]
    max_times = get_intersection_times(integrators)
    hull = Hull(integrators, max_times)
    density_grid = update_density_grid(wi.density_grid, wi.update_grid_flag, integrators, max_times, hull)
    velocity_grid = update_velocity_grid(wi.velocity_grid, wi.update_grid_flag, integrators, max_times, hull)
    return WindInterpolator(
        hull,
        density_grid,
        velocity_grid,
        wi.vacuum_density,
        wi.n_timesteps,
        wi.update_grid_flag,
    )
end

function update_wind_interpolator(wi::WindInterpolator, dgrid::DensityGrid)
    return WindInterpolator(
        wi.wind_hull,
        dgrid,
        wi.velocity_grid,
        wi.vacuum_density,
        wi.n_timesteps,
        wi.update_grid_flag
    )
end

function update_wind_interpolator(wi::WindInterpolator, vgrid::VelocityGrid)
    return WindInterpolator(
        wi.wind_hull,
        wi.density_grid,
        vgrid,
        wi.vacuum_density,
        wi.n_timesteps,
        wi.update_grid_flag
    )
end

function get_density(wi::WindInterpolator, r, z)
    if is_point_in_wind(wi, [r, z])
        return wi.grid.interpolator(r, z)
    else
        return wi.vacuum_density
    end
end
get_density(wi::WindInterpolator, point) = get_density(wi, point[1], point[2])

function GridIterator(interpolator::WindInterpolator, ri, zi, rf, zf)
    return GridIterator(
        interpolator.density_grid.r_range,
        interpolator.density_grid.z_range,
        ri,
        zi,
        rf,
        zf,
    )
end
