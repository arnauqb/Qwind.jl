using PyCall
import ConcaveHull, Interpolations, Sundials
export WindInterpolator, get_density, update_interpolator, scipy_interpolate

const scipy_interpolate = PyNULL()

function __init__()
    copy!(scipy_interpolate, pyimport_conda("scipy.interpolate", "scipy"))
end


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
        density_grid = construct_density_grid(nr, nz)
        velocity_grid = construct_velocity_grid(nr, nz)
    else
        r0 = [integ.p.r0 for integ in integrators]
        @info "Constructing wind hull..."
        hull = construct_wind_hull(integrators)
        flush()
        @info "Constructing interpolation grid..."
        flush()
        r, z, vr, vz, n = reduce_integrators(integrators, n_timesteps = 1000)
        density_grid = construct_density_grid(r, z, n, r0, hull, nr = nr, nz = nz)
        velocity_grid = construct_velocity_grid(r, z, vr, vz, r0, hull, nr = nr, nz = nz)
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
    @info "Constructing wind hull..."
    flush()
    new_hull = construct_wind_hull(integrators)
    @info "Constructing interpolation grid..."
    r0s = [integ.p.r0 for integ in integrators]
    r, z, vr, vz, n = reduce_integrators(integrators, n_timesteps = 1000)
    flush()
    density_grid = update_density_grid(wi.density_grid, new_hull, r, z, n, r0s)
    velocity_grid = construct_velocity_grid(
        r,
        z,
        vr,
        vz,
        r0s,
        new_hull,
        nr = wi.velocity_grid.nr,
        nz = wi.velocity_grid.nz,
    )
    return WindInterpolator(
        new_hull,
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
