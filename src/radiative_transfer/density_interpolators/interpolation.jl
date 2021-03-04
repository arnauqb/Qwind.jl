#__precompile__()
using PyCall, BasisFunctionExpansions
import ConcaveHull, Interpolations
export WindInterpolator, get_density

#scipy_interpolate = PyNULL()

#function __init__()
#    copy!(scipy_interpolate, pyimport_conda("scipy.interpolate", "scipy"))
#end
#
struct WindInterpolator{T} <: GridInterpolator{T}
    grid::InterpolationGrid{T}
    hull::Union{ConcaveHull.Hull,Nothing}
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
        grid = construct_interpolation_grid(nr, nz)
    else
        r0 = [integ.p.r0 for integ in integrators]
        @info "Constructing wind hull..."
        hull = construct_wind_hull(integrators)
        flush()
        @info "Constructing interpolation grid..."
        flush()
        r, z, n = reduce_integrators(integrators, n_timesteps = 1000)
        grid = construct_interpolation_grid(r, z, n, r0, hull, nr = nr, nz = nz)
    end
    return WindInterpolator(grid, hull, vacuum_density, n_timesteps)
end

function construct_wind_hull(r::Vector{Float64}, z::Vector{Float64}, r0::Vector{Float64})
    r_log = log10.(r)
    z_log = log10.(z)
    r0_log = log10.(r0)
    r_range_log = log10.(range(minimum(r), maximum(r), step = 5))
    r_range_log = sort(vcat(r_range_log, r0_log))
    z_range_log = range(max(minimum(z_log), -8), maximum(z_log), length = 1000)

    zmax_hull = -Inf .* ones(length(r_range_log))
    rmax_hull = -Inf .* ones(length(z_range_log))
    rmin_hull = Inf .* ones(length(z_range_log))

    # top part of the wind
    for (rp_log, zp_log) in zip(r_log, z_log)
        ridx = searchsorted_nearest(r_range_log, rp_log)
        zidx = searchsorted_nearest(z_range_log, zp_log)

        # get highest position of the wind
        zmax_hull[ridx] = max(zmax_hull[ridx], zp_log)

        # get furthest out part of the wind
        rmax_hull[zidx] = max(rmax_hull[zidx], rp_log)

        # get innermost out part of the wind
        rmin_hull[zidx] = min(rmin_hull[zidx], rp_log)
    end

    r_hull = r0_log
    z_hull = -10 .* ones(length(r_hull))

    # filter unassigned
    mask = zmax_hull .!= -Inf
    r_hull = vcat(r_hull, r_range_log[mask])
    z_hull = vcat(z_hull, zmax_hull[mask])

    mask = rmax_hull .!= -Inf
    r_hull = vcat(r_hull, rmax_hull[mask])
    z_hull = vcat(z_hull, z_range_log[mask])

    mask = rmin_hull .!= Inf
    r_hull = vcat(r_hull, rmin_hull[mask])
    z_hull = vcat(z_hull, z_range_log[mask])

    points = []
    for (rp, zp) in zip(r_hull, z_hull)
        push!(points, [rp, zp])
    end
    points = reduce(hcat, points)
    points = round.(points, sigdigits = 3)
    points = unique(points, dims = 2)
    points = [[points[1, i], points[2, i]] for i = 1:size(points)[2]]
    #return points
    hull = ConcaveHull.concave_hull(points)
    return hull#, points
end

function get_dense_line_positions(integrators; n_timesteps = 1000)
    r = Float64[]
    z = Float64[]
    for integ in integrators
        rp, zp, _, _, _, _ = get_dense_solution_from_integrator(integ, n_timesteps)
        r = vcat(r, rp)
        z = vcat(z, zp)
    end
    return r, z
end

function construct_wind_hull(integrators)
    r0 = [integ.p.r0 for integ in integrators]
    hull = nothing
    for nt in [100, 1000, 5000]
        @info "Trying wind hull with $nt time steps..."
        flush()
        r, z = get_dense_line_positions(integrators, n_timesteps = nt)
        hull = construct_wind_hull(r, z, r0)
        hull.converged && break
    end
    if hull === nothing 
        error("Cannot construct wind hull!")
    end
    return hull
end

function is_point_in_wind(hull::ConcaveHull.Hull, point)
    return ConcaveHull.in_hull(log10.(point), hull)
end
is_point_in_wind(wi::WindInterpolator, point) = is_point_in_wind(wi.hull, point)
is_point_in_wind(hull::ConcaveHull.Hull, r, z) = is_point_in_wind(hull, [r, z])

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


function reduce_integrators(integrators; n_timesteps = 10000)
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

function construct_interpolation_grid(nr, nz)
    r_range = zeros(2)
    z_range = zeros(2)
    density_grid = [[1e2, 1e2] [1e2, 1e2]]
    return InterpolationGrid(r_range, z_range, density_grid, nr, nz)
end

function get_density_interpolator(
    r::Vector{Float64},
    z::Vector{Float64},
    n::Vector{Float64},
)
    r_log = log10.(r)
    z_log = log10.(z)
    log_n = log10.(n)
    points = hcat(r_log, z_log)
    scipy_interpolate = pyimport("scipy.interpolate")
    linear_int = scipy_interpolate.LinearNDInterpolator(points, log_n, fill_value = 2)
    return linear_int
end

function construct_interpolation_grid(
    r::Vector{Float64},
    z::Vector{Float64},
    n::Vector{Float64},
    r0s::Vector{Float64},
    hull;
    nr = "auto",
    nz = 50,
)
    interpolator = get_density_interpolator(r, z, n)
    r_range, z_range = get_spatial_grid(r, z, r0s, nr, nz)
    r_range_log = log10.(r_range)
    z_range_log = log10.(z_range)
    #r_range_grid_log = r_range_log .* ones(length(z_range_log))'
    #z_range_grid_log = z_range_log' .* ones(length(r_range_log))
    #density_grid = 10 .^ linear_int(log10.(r_range_grid), log10.(z_range_grid));
    #density_grid = @showprogress pmap(
    #    z_log -> 10 .^ linear_int(r_range_log, z_log),
    #    z_range_log,
    #    batch_size = Int(round(length(z_range) / nprocs())),
    #)
    #density_grid = reduce(hcat, density_grid)
    #density_grid =
    #    10 .^ scipy_interpolate.griddata(
    #        (r_log, z_log),
    #        log_n,
    #        (r_range_grid_log, z_range_grid_log),
    #        method = "linear",
    #        fill_value = 2,
    #    )
    # remove points outside the wind
    density_grid = 1e2 .* ones((length(r_range), length(z_range)))
    for (i, r) in enumerate(r_range)
        for (j, z) in enumerate(z_range)
            point = [r, z]
            if !is_point_in_wind(hull, point)
                density_grid[i, j] = 1e2
            else
                density_grid[i, j] = 10 .^ interpolator(log10(r), log10(z))[1]
            end
        end
    end
    # add z = 0 line
    density_grid = [density_grid[:, 1] density_grid]
    pushfirst!(z_range, 0.0)
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
    #z_range = pushfirst!(z_range, 0.0)
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
    if maximum(wi.grid.grid) == wi.vacuum_density
        # first iteration, do not average
        return WindInterpolator(
            integrators,
            nr = wi.grid.nr,
            nz = wi.grid.nz,
            vacuum_density = wi.vacuum_density,
            n_timesteps = wi.n_timesteps,
        )
    end
    old_grid = wi.grid
    @info "Constructing wind hull..."
    new_hull = construct_wind_hull(integrators)
    @info "Constructing interpolation grid..."
    r0s = [integ.p.r0 for integ in integrators]
    r, z, n = reduce_integrators(integrators, n_timesteps = 1000)
    r_range, z_range = get_spatial_grid(r, z, r0s, wi.grid.nr, wi.grid.nz)
    r_range_log = log10.(r_range)
    z_range_log = log10.(z_range)
    interpolator = get_density_interpolator(r, z, n)
    density_grid = 1e2 .* ones((length(r_range), length(z_range)))
    for (i, r) in enumerate(r_range)
        for (j, z) in enumerate(z_range)
            point = [r, z]
            if !is_point_in_wind(new_hull, point)
                density_grid[i, j] = 1e2
            else
                density_grid[i, j] =
                    10 .^ (
                        (
                            interpolator(log10(r), log10(z))[1] +
                            log10(get_density(old_grid, r, z))
                        ) / 2.0
                    )
            end
        end
    end
    # add z = 0 line
    density_grid = [density_grid[:, 1] density_grid]
    pushfirst!(z_range, 0.0)
    grid = InterpolationGrid(r_range, z_range, density_grid, wi.grid.nr, wi.grid.nz)
    return WindInterpolator(grid, new_hull, wi.vacuum_density, wi.n_timesteps)
end

function get_density(wi::WindInterpolator, iterator::GridIterator)
    return wi.grid.grid[iterator.current_r_idx, iterator.current_z_idx]
end
