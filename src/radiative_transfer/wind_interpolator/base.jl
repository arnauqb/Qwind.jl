using PyCall
import ConcaveHull, Interpolations, Sundials
export WindInterpolator, get_density, update_interpolator

struct WindInterpolator{T} <: Interpolator{T}
    wind_hull::Union{ConcaveHull.Hull,Nothing}
    density_grid::InterpolationGrid{T}
    velocity_grid::InterpolationGrid{T}
    vacuum_density::T
    n_timesteps::Int
end

function WindInterpolator(
    integrators;
    nr = "auto",
    nz = 50,
    vacuum_density = 1e2,
    n_timesteps = 1000,
)
    if nr != "auto"
        nr = Int(nr)
    end
    nz = Int(nz)
    if integrators === nothing
        hull = nothing
        density_grid = DensityGrid(nr, nz, vacuum_density)
        velocity_grid = VelocityGrid(nr, nz, 0.0)
    else
        r0 = [integ.p.r0 for integ in integrators]
        max_times = get_intersection_times(integrators)
        hull = Hull(integrators, max_times)
        density_grid = DensityGrid(integrators, max_times, hull, nr=nr, nz=nz)
        velocity_grid = VelocityGrid(integrators, max_times, hull, nr=nr, nz=nz)
    end
    return WindInterpolator(hull, density_grid, velocity_grid, vacuum_density, n_timesteps)
end

function update_interpolator(wi::WindInterpolator, integrators)
    if maximum(wi.density_grid.grid) == wi.vacuum_density
        # first iteration, do not average
        return WindInterpolator(
            integrators,
            nr = wi.density_grid.nr,
            nz = wi.density_grid.nz,
            vacuum_density = wi.vacuum_density,
            n_timesteps = wi.n_timesteps,
        )
    end
    r0 = [integ.p.r0 for integ in integrators]
    max_times = get_intersection_times(integrators)
    hull = Hull(integrators, max_times)
    density_grid = update_density_grid(wi.density_grid, integrators, max_times, hull)
    velocity_grid = update_velocity_grid(wi.velocity_grid, integrators, max_times, hull)
    return WindInterpolator(
        hull,
        density_grid,
        velocity_grid,
        wi.vacuum_density,
        wi.n_timesteps,
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
