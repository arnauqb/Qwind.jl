using Dierckx
using PyCall
export LogGrid, update_density_interpolator
scipy_interp = pyimport("scipy.interpolate")

struct LogGrid <: GridInterpolator
    grid::InterpolationGrid
    kdtree::Union{KDTree,Nothing}
    vacuum_density::Float64
    kx::Int
    ky::Int
    s::Float64
    n_timesteps::Int
    function LogGrid(
        r_range,
        z_range,
        grid;
        kdtree,
        vacuum_density = 1e2,
        kx = 3,
        ky = 3,
        s = 0,
        n_timesteps = 10000,
    )
        return new(
            InterpolationGrid(r_range, z_range, grid),
            kdtree,
            vacuum_density,
            kx,
            ky,
            s,
            n_timesteps,
        )
    end
end

function LogGrid(
    integrators;
    nr = 250,
    nz = 250,
    kx = 3,
    ky = 3,
    s = 0.0,
    vacuum_density = 1e2,
    n_timesteps = 10000,
)
    if integrators === nothing
        kdtree = nothing
    else
        kdtree = create_wind_kdtree(integrators, n_timesteps, vacuum_density)
    end
    return LogGrid(
        kdtree,
        nr,
        nz,
        kx = kx,
        ky = ky,
        s = s,
        vacuum_density = vacuum_density,
        n_timesteps = n_timesteps,
    )
end

function LogGrid(
    kdtree::Union{KDTree,Nothing},
    nr,
    nz;
    kx = 3,
    ky = 3,
    s = 0.0,
    vacuum_density = 1e2,
    n_timesteps = 10000,
)
    if kdtree === nothing
        r_range = zeros(nr)
        z_range = zeros(nz)
        density_grid = nothing 
    else
        r_range, z_range = get_spatial_grid(kdtree, nr, nz)
        #z_range = z_range[1:end-1]
        grid_r = r_range' .* ones(length(z_range));
        grid_z = z_range' .* ones(length(r_range));
        #density_grid = 10 .^ scipy_interp.griddata(
        #    (kdtree.r, kdtree.z),
        #    log10.(kdtree.n),
        #    (grid_r, grid_z),
        #    fill_value = 2,
        #    method="nearest"
        #)
        #spline = Spline2D(kdtree.r, kdtree.z, log10.(kdtree.n))
        #density_grid =
        #    10 .^ reshape(
        #        [get_density(spline, kdtree, r, z) for z in z_range for r in r_range],
        #        (nr, nz),
        #    )
        # rbf
        #rbfi = scipy_interp.Rbf(kdtree.r, kdtree.z, kdtree.n)
        #density_grid = 10 .^ reshape(rbfi(grid_r, grid_z), (nr, nz))
    end
    return LogGrid(
        r_range,
        z_range,
        density_grid,
        kdtree = kdtree,
        vacuum_density = vacuum_density,
        kx = kx,
        ky = ky,
        s = s,
        n_timesteps = n_timesteps,
    )
end

function get_density(spline::Spline2D, kdtree::KDTree, r, z)
    rc, zc, zmax, z0, width, density, distance = get_closest_point(kdtree, r, z)
    if is_point_outside_wind(z, zmax, z0, distance, width)
        return kdtree.vacuum_density
    else
        return Dierckx.evaluate(spline, r, z)
    end
end

function update_density_interpolator(interpolator::LogGrid, integrators)
    new_interpolator = LogGrid(
        integrators,
        nr = interpolator.grid.nr,
        nz = interpolator.grid.nz,
        kx = interpolator.kx,
        ky = interpolator.ky,
        s = interpolator.s,
        vacuum_density = interpolator.vacuum_density,
        n_timesteps = interpolator.n_timesteps,
    )
    return new_interpolator
end
