#__precompile__()
using PyCall
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
    n_timesteps = 10000,
)
    if nr != "auto"
        nr = Int(nr)
    end
    nz = Int(nz)
    if integrators === nothing
        hull = nothing
        grid = construct_interpolation_grid(nr, nz)
    else
        r, z, n = reduce_integrators(integrators, n_timesteps=n_timesteps)
        r0 = [integ.p.r0 for integ in integrators]
        hull = construct_wind_hull(r, z, r0)
        grid = construct_interpolation_grid(r, z, n, r0, hull, nr = nr, nz = nz)
    end
    return WindInterpolator(grid, hull, vacuum_density, n_timesteps)
end

function construct_wind_hull(r::Vector{Float64}, z::Vector{Float64}, r0::Vector{Float64})
    ri = minimum(r)
    rf = maximum(r)
    current_r = ri
    r_range = [ri]
    while current_r < rf
        r0_idx = searchsorted_first(r0, current_r) + 1
        if r0_idx > length(r0)
            deltar = 10
        else
            deltar0 = r0[r0_idx] - current_r
            deltar = min(1, deltar0)
        end
        current_r += deltar
        @assert deltar > 0
        push!(r_range, current_r)
    end
    zmax_hull = -1 .* ones(length(r_range))
    for (rp, zp) in zip(r, z)
        r_idx = searchsorted_nearest(r_range, rp)
        zmax_hull[r_idx] = max(zmax_hull[r_idx], zp)
    end
    mask = zmax_hull .>= 0
    zmax_hull = zmax_hull[mask]
    r_range = r_range[mask]
    zmin_hull = zeros(length(r0))
    z_hull = vcat(zmax_hull, zmin_hull)
    r_hull = vcat(r_range, r0)
    points = []
    for (rp, zp) in zip(r_hull, z_hull)
        push!(points, [rp, zp])
    end
    hull = ConcaveHull.concave_hull(points)
    return hull
end

function get_dense_line_positions(integrators; n_timesteps=1000)
    r = Float64[]
    z = Float64[]
    for integ in integrators
        rp, zp, _, _, _, _ = get_dense_solution_from_integrator(integ, n_timesteps)
        r = vcat(r, rp)
        z = vcat(z, zp)
    end
    return r, z
end

function construct_wind_hull(integrators; n_timesteps = 1000)
    r, z = get_dense_line_positions(integrators, n_timesteps=n_timesteps)
    r0 = [integ.p.r0 for integ in integrators]
    return construct_wind_hull(r, z, r0)
end

function is_point_in_wind(hull::ConcaveHull.Hull, point)
    return ConcaveHull.in_hull(point, hull)
end
is_point_in_wind(wi::WindInterpolator, point) = is_point_in_wind(wi.hull, point)
is_point_in_wind(hull::ConcaveHull.Hull, r, z) = is_point_in_wind(hull, [r,z])

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

function construct_interpolation_grid(nr, nz)
    r_range = zeros(2)
    z_range = zeros(2)
    density_grid = [[1e2, 1e2] [1e2, 1e2]]
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
