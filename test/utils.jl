export create_test_quadtree

constant_density(r, z) = 1e8 * ones(length(r))
linear_density(r, z) = 1e8 .* (r + 2 * z)
powerlaw_density_1(r, z) = 1e8 ./ (r + z)
powerlaw_density_2(r, z) = 1e8 .* (r + z)
exponential_density(r, z) = 1e8 .* exp.(-(r + z) / 20)
density_functions = [
    constant_density,
    linear_density,
    powerlaw_density_1,
    powerlaw_density_2,
    exponential_density,
]

function create_test_quadtree(
    density_function;
    r_range,
    z_range,
    atol = 5e-4,
    rtol = 1e-3,
    cell_min_size = 1e-6,
    width_range=nothing,
)
    bh = BlackHole(1e8 * M_SUN, 0.5, 0)
    zmax = maximum(z_range) * ones(length(r_range))
    if width_range === nothing
        width_range = 300 .* ones(length(zmax))
    end
    density = density_function(r_range, z_range)
    kdtree = create_wind_kdtree(r_range, z_range, zmax, width_range, density)
    quadtree = create_and_refine_quadtree(
        kdtree,
        Rg = bh.Rg,
        atol = atol,
        rtol = rtol,
        cell_min_size = cell_min_size,
    )
    return quadtree, kdtree, bh
end
