using ImageFiltering
using Qwind

export SmoothedGrid,
    get_density, get_density_points, r_min, r_max, z_min, z_max, update_density_interpolator

struct SmoothedGrid <: GridInterpolator
    grid::InterpolationGrid
    kdtree::Union{KDTree,Nothing}
    vacuum_density::Float64
    kernel_size::Int
    n_timesteps::Int
    function SmoothedGrid(
        r_range,
        z_range,
        grid;
        kdtree,
        vacuum_density = 1e2,
        kernel_size = 1,
        n_timesteps = 10000,
    )
        return new(
            InterpolationGrid(r_range, z_range, grid),
            kdtree,
            vacuum_density,
            kernel_size,
            n_timesteps,
        )
    end
end

function SmoothedGrid(
    integrators;
    nr = 250,
    nz = 250,
    kernel_size = 1,
    vacuum_density = 1e2,
    n_timesteps = 10000,
)
    if integrators === nothing
        kdtree = nothing
    else
        kdtree = create_wind_kdtree(integrators, n_timesteps, vacuum_density)
    end
    return SmoothedGrid(
        kdtree,
        nr,
        nz,
        kernel_size = kernel_size,
        vacuum_density = vacuum_density,
        n_timesteps = n_timesteps,
    )
end

function SmoothedGrid(
    kdtree::Union{KDTree,Nothing},
    nr,
    nz;
    kernel_size = 10,
    vacuum_density = 1e2,
    n_timesteps = 10000,
)
    if kdtree === nothing
        r_range = zeros(nr) #[0.0, 1.0]
        z_range = zeros(nz) #[0.0, 1.0]
        density_grid = nothing #vacuum_density .* ones((2, 2))
    else
        r_range, z_range = get_spatial_grid(kdtree, nr, nz)
        points = hcat([[r, z] for r in r_range for z in z_range]...)
        density_grid = reshape(get_density_points(kdtree, points), (nz, nr))'
        density_grid = imfilter(density_grid, Kernel.gaussian(kernel_size))
    end
    return SmoothedGrid(
        r_range,
        z_range,
        density_grid,
        kdtree=kdtree,
        vacuum_density=vacuum_density,
        kernel_size=kernel_size,
        n_timesteps=n_timesteps,
    )
end

function update_density_interpolator(interpolator::SmoothedGrid, integrators)
    new_interpolator = SmoothedGrid(
        integrators,
        nr = interpolator.grid.nr,
        nz = interpolator.grid.nz,
        kernel_size = interpolator.kernel_size,
        vacuum_density = interpolator.vacuum_density,
        n_timesteps = interpolator.n_timesteps,
    )
    return new_interpolator
end

