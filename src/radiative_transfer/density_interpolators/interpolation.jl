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
    grid = construct_interpolation_grid(
        integrators,
        hull,
        n_timesteps = n_timesteps,
        nr = nr,
        nz = nz,
    )
    return WindInterpolator(grid, hull, vacuum_density, n_timesteps)
end

function construct_wind_hull(r::Vector{Float64}, z::Vector{Float64})
    r_hull = range(minimum(r), maximum(r), step = 1)
    zmax_hull = zeros(length(r_hull))
    zmin_hull = zeros(length(r_hull))
    for (rp, zp) in zip(r, z)
        r_idx = searchsorted_nearest(r_hull, rp)
        zmax_hull[r_idx] = max(zmax_hull[r_idx], zp)
        zmin_hull[r_idx] = min(zmin_hull[r_idx], zp)
    end
    z_hull = vcat(zmax_hull, zmin_hull)
    points = []
    for (rp, zp) in zip(r_hull, z_hull)
        #zp == 0 && continue
        push!(points, [rp, zp])
    end
    hull = ConcaveHull.concave_hull(points)
    return hull
end

function construct_wind_hull(integrators; n_timesteps = 10000)
    r = Float64[]
    z = Float64[]
    for integ in integrators
        rp, zp, _, _, _, _ = get_dense_solution_from_integrator(integ, n_timesteps)
        r = vcat(r, rp)
        z = vcat(z, zp)
    end
    return construct_wind_hull(r, z)
end

function is_point_in_wind(hull::ConcaveHull.Hull, point)
    return ConcaveHull.in_hull(point, hull)
end
is_point_in_wind(wi::WindInterpolator, point) = is_point_in_wind(wi.hull, point)

function reduce_line(r::Vector{<:Number}, z::Vector{<:Number}, n::Vector{<:Number})
    rs = [r[1]]
    zs = [z[1]]
    ns = [n[1]]
    for (rp, zp, np) in zip(r, z, n)
        if (zp > zs[end])
            push!(rs, rp)
            push!(zs, zp)
            push!(ns, np)
        end
    end
    return rs, zs, ns
end

function reduce_line(integrator::Sundials.IDAIntegrator; n_timesteps = 10000)
    r, z, zmax, z0, width, n = get_dense_solution_from_integrator(integrator, n_timesteps)
    return reduce_line(r, z, n)
end


function reduce_integrators(integrators; n_timesteps)
    rs = Float64[]
    zs = Float64[]
    ns = Float64[]
    for integ in integrators
        r, z, n = reduce_line(integ, n_timesteps = n_timesteps)
        rs = vcat(rs, r)
        zs = vcat(zs, z)
        ns = vcat(ns, n)
    end
    return rs, zs, ns
end

function construct_interpolation_grid(
    integrators::Union{Vector{<:Sundials.IDAIntegrator}, Nothing},
    hull::Union{ConcaveHull.Hull, Nothing};
    n_timesteps,
    nr,
    nz,
)
    if integrators === nothing
        return construct_interpolation_grid(nr, nz)
    else
        rs, zs, ns = reduce_integrators(integrators, n_timesteps = n_timesteps)
        r0s = [integrator.p.r0 for integrator in integrators]
        return construct_interpolation_grid(rs, zs, ns, r0s, hull, nr = nr, nz = nz)
    end
end

function construct_interpolation_grid(nr, nz)
    r_range = zeros(2)
    z_range = zeros(2)
    density_grid = nothing
    return InterpolationGrid(r_range, z_range, density_grid, nr, nz)
end

function construct_interpolation_grid(
    r::Vector{Float64},
    z::Vector{Float64},
    n::Vector{Float64},
    r0s::Vector{Float64},
    hull;
    nr="auto",
    nz=50,
)
    log_n = log10.(n)
    r_range, z_range = get_spatial_grid(r, z, r0s, nr, nz)
    #r_range_grid = r_range .* ones(length(z_range))'
    #z_range_grid = z_range' .* ones(length(r_range))
    scipy_interpolate = pyimport("scipy.interpolate")
    points = hcat(r, z)
    linear_int = scipy_interpolate.LinearNDInterpolator(points, log_n, fill_value=2)
    #density_grid = zeros((length(r_range), length(z_range)))
    #@showprogress for (i, r) in enumerate(r_range)
    #    for (j, z) in enumerate(z_range)
    #        density_grid[i,j] = 10.0 ^ linear_int(r, z)[1]
    #    end
    #end
    density_grid = @showprogress pmap(
        z -> 10 .^ linear_int(r_range, z),
        z_range,
        batch_size = Int(round(length(z_range) / nprocs())),
    )
    density_grid = hcat(density_grid...)
    #density_grid =
    #    10 .^ scipy_interpolate.griddata(
    #        (r, z),
    #        log_n,
    #        (r_range_grid, z_range_grid),
    #        method = "linear",
    #        fill_value = 2,
    #    )
    # remove points outside the wind
    for (i, r) in enumerate(r_range)
        for (j, z) in enumerate(z_range)
            point = [r, z]
            if !is_point_in_wind(hull, point)
                density_grid[i, j] = 1e2
            end
        end
    end
    grid = InterpolationGrid(r_range, z_range, density_grid, nr, nz)
    return grid
end

function get_spatial_grid(
    r::Vector{Float64},
    z::Vector{Float64},
    r0s::Vector{Float64} = nothing,
    nr = "auto",
    nz = 50,
)
    r_min = max(6.0, minimum(r))
    z_min = max(minimum(z), 1e-6)
    r_max = min(1e4, maximum(r))
    z_max = min(1e4, maximum(z))
    if nr == "auto" 
        r_range = r0s
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
    if is_point_in_wind(wi, [r, z])
        return wi.grid.interpolator(r, z)
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

function get_density(wi::WindInterpolator, iterator::GridIterator)
    return wi.grid.grid[iterator.current_r_idx, iterator.current_z_idx]
end
