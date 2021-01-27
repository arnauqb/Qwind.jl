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

function create_test_bh()
    bh = BlackHole(1e8 * M_SUN, 0.5, 0)
    bh
end

function create_test_kdtree(density_function; r_range, z_range, width_range)
    zmax = maximum(z_range) * ones(length(r_range))
    if width_range === nothing
        width_range = 300 .* ones(length(zmax))
    end
    density = density_function(r_range, z_range)
    kdtree = create_wind_kdtree(r_range, z_range, zmax, width_range, density)
    kdtree
end

function create_test_quadtree(
    density_function;
    r_range,
    z_range,
    atol = 5e-4,
    rtol = 1e-3,
    cell_min_size = 1e-6,
    width_range = nothing,
)
    kdtree = create_test_kdtree(
        density_function,
        r_range = r_range,
        z_range = z_range,
        width_range = width_range,
    )
    bh = create_test_bh()
    quadtree = create_and_refine_quadtree(
        kdtree,
        Rg = bh.Rg,
        atol = atol,
        rtol = rtol,
        cell_min_size = cell_min_size,
    )
    return quadtree
end

unzip(a) = map(x->getfield.(a, x), fieldnames(eltype(a)))

function remove_close_elements(args...; digits=4)
    ret = [array[1] for array in args]
    for i in 2:length(args[1])
        for (j, array) in enumerate(args)
            element = trunc(array[i], digits=digits)
            if ret[i-1, j] == element 
                continue
            else
                ret[i, j] = element
            end
        end
    end
end
