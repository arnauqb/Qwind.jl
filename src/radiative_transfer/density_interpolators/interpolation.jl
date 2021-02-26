#__precompile__()
using PyCall
import ConcaveHull, Interpolations
export WindInterpolator, get_density

#scipy_interpolate = PyNULL() 

#function __init__()
#    copy!(scipy_interpolate, pyimport_conda("scipy.interpolate", "scipy"))
#end
#
struct WindInterpolator <: GridInterpolator
    grid::InterpolationGrid
    hull::Union{ConcaveHull.Hull,Nothing}
    vacuum_density::Float64
    n_timesteps::Int
end

function WindInterpolator(
    integrators;
    nr = "auto",
    nz = 50,
    vacuum_density = 1e2,
    n_timesteps = 10000,
)
    if nr != "auto"
        nr = Int(nr)
    end
    nz = Int(nz)
    if integrators === nothing
        hull = nothing
    else
        hull = construct_wind_hull(integrators, n_timesteps = n_timesteps)
    end
    grid = construct_interpolation_grid(integrators, hull, n_timesteps = n_timesteps, nr=nr, nz=nz)
    return WindInterpolator(
        grid,
        hull,
        vacuum_density,
        n_timesteps,
    )
end

function construct_wind_hull(integrators; n_timesteps = 10000)
    r_min = minimum([minimum(integ.p.data[:r]) for integ in integrators])
    r_max = maximum([maximum(integ.p.data[:r]) for integ in integrators])
    r_range = range(r_min, r_max, step = 1)
    z_range = zeros(length(r_range))
    for integ in integrators
        r, z, zmax, z0, width, n = get_dense_solution_from_integrator(integ, n_timesteps)
        for (rp, zp) in zip(r, z)
            r_idx = searchsorted_nearest(r_range, rp)
            z_range[r_idx] = max(z_range[r_idx], zp)
        end
    end
    points = [[integ.p.r0, 0] for integ in integrators]
    for (rp, zp) in zip(r_range, z_range)
        zp == 0 && continue
        push!(points, [rp, zp])
    end
    hull = ConcaveHull.concave_hull(points)
    return hull
end

function is_point_in_wind(hull::ConcaveHull.Hull, point)
    return ConcaveHull.in_hull(point, hull)
end
is_point_in_wind(wi::WindInterpolator, point) = is_point_in_wind(wi.hull, point)

function reduce_line(integrator; n_timesteps = 10000)
    r, z, zmax, z0, width, n = get_dense_solution_from_integrator(integrator, n_timesteps)
    rs = [r[1]]
    zs = [z[1]]
    ns = [n[1]]
    for (rp, zp, np) in zip(r, z, n)
        #z_ratio = abs(zp / zs[end])
        if (zp > zs[end])# && (z_ratio < 0.95 || z_ratio > 1.05)
            push!(rs, rp)
            push!(zs, zp)
            push!(ns, np)
        end
    end
    return rs, zs, ns
end

function reduce_integrators(integrators; n_timesteps)
    rs = []
    zs = []
    ns = []
    for integ in integrators
        r, z, n = reduce_line(integ, n_timesteps = n_timesteps)
        rs = vcat(rs, r)
        zs = vcat(zs, z)
        ns = vcat(ns, n)
    end
    return rs, zs, ns
end

function construct_interpolation_grid(integrators, hull; n_timesteps, nr, nz)
    if integrators === nothing
        r_range = zeros(2)
        z_range = zeros(2)
        density_grid = nothing
        interpolator = (r,z) -> 1e2
    else
        rs, zs, ns = reduce_integrators(integrators, n_timesteps=n_timesteps)
        ns = log10.(ns)
        r_range, z_range = get_spatial_grid(integrators, nr, nz)
        r_range_grid = r_range .* ones(length(z_range))'
        z_range_grid = z_range' .* ones(length(r_range))
        scipy_interpolate = pyimport("scipy.interpolate")
        density_grid =
            10 .^ scipy_interpolate.griddata(
                (rs, zs),
                ns,
                (r_range_grid, z_range_grid),
                method = "linear",
                fill_value = 2,
            )
        # remove points outside the wind
        for (i, r) in enumerate(r_range)
            for (j,z) in enumerate(z_range)
                point = [r,z]
                if !is_point_in_wind(hull, point)
                    density_grid[i,j] = 1e2
                end
            end
        end
        interpolator = Interpolations.interpolate((r_range, z_range), density_grid, Gridded(Linear()))
        interpolator = Interpolations.extrapolate(interpolator, 1e2)
    end
    grid = InterpolationGrid(r_range, z_range, density_grid, interpolator, nr, nz)
    return grid
end

function get_spatial_grid(integrators, nr, nz)
    r_min = Inf
    r_max = 0
    z_min = Inf
    z_max = 0
    for integrator in integrators
        r_min = min(r_min, integrator.p.r0)
        r_max = max(r_max, min(maximum(integrator.p.data[:r] * (1 + integrator.p.lwnorm)), 5e3))
        z_min = min(z_min, integrator.p.z0)
        z_max = max(z_max, maximum(integrator.p.data[:z]))
    end
    r_min = max(6, r_min)
    z_min = max(1e-6, z_min)
    if nr == "auto"
        r_range = [integ.p.r0 for integ in integrators]
        additional_r = collect(range(r_range[end], r_max, step = 100)[2:end])
        r_range = vcat(r_range, additional_r)
    else
        r_range = 10 .^ range(log10(r_min), log10(r_max), length = nr)
    end
    z_range = 10 .^ range(log10(z_min), log10(z_max), length = nz - 1)
    z_range = pushfirst!(z_range, 0.0)
    r_range = round.(r_range, digits = 7)
    z_range = round.(z_range, digits = 7)
    points = hcat([[r, z] for r in r_range for z in z_range]...)
    return r_range, z_range
end

function get_density(wi::WindInterpolator, r, z)
    if is_point_in_wind(wi, [r,z])
        return wi.grid.interpolator(r,z)
    else
        return wi.vacuum_density
    end
end
get_density(wi::WindInterpolator, point) = get_density(wi, point[1], point[2])

function update_density_interpolator(wi::WindInterpolator, integrators)
    return WindInterpolator(
        integrators;
        nr = wi.grid.nr,
        nz = wi.grid.nz,
        vacuum_density = wi.vacuum_density,
        n_timesteps = wi.n_timesteps,
    )
end

function get_density(
    wi::WindInterpolator,
    iterator::GridIterator,
)
    return wi.grid.grid[iterator.current_r_idx, iterator.current_z_idx]
end
