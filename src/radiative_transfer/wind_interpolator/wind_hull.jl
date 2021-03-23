import ConcaveHull
import ConcaveHull: Hull
export DenseIntegrator, DenseIntegrators, Hull

function Hull(r::Vector{Float64}, z::Vector{Float64}, r0::Vector{Float64}; sigdigits=6)
    # this stores the maximum heights of each streamline, roughly
    z_max = zeros(length(r0))
    for (rp, zp) in zip(r, z)
        if rp > r0[end] + 1
            continue
        end
        ridx = searchsorted_nearest(r0, rp)
        z_max[ridx] = max(z_max[ridx], zp)
    end
    ridx = searchsorted_nearest(r0, r[end])
    r0 = r0[z_max .> 0]
    z_max = z_max[z_max .> 0]

    # then get points in between trajectories, this avoids
    # the hull being too tight between trajectories
    r_range = range(minimum(r0), maximum(r0), step=1)
    z_range = zeros(length(r_range))
    for (i, rp) in enumerate(r_range)
        ridx = searchsorted_nearest(r0, rp)
        z_range[i] = max(z_range[i], z_max[ridx])
    end
    r = vcat(r, r_range)
    z = vcat(z, z_range)

    # remove points that are too close to each other
    points = hcat(r, z)
    points = round.(points, sigdigits = sigdigits)
    points = unique(points, dims = 1)
    points = [[points[i, 1], points[i, 2]] for i = 1:size(points)[1]]
    @info "Constructing wind hull..."
    flush()
    hull = ConcaveHull.concave_hull(points)
    if !hull.converged
        error("Hull did not converge!")
    end
    @info "Done"
    return hull#, points
end

function Hull(integrators::Vector{<:Sundials.IDAIntegrator}, max_times; hull_sigdigits=6)
    r0 = [integ.p.r0 for integ in integrators]
    integrators_interpolated_linear = interpolate_integrators(
        integrators,
        max_times = max_times,
        n_timesteps = 50,
        log = false,
    )
    r, z, _, _, _ = reduce_integrators(integrators_interpolated_linear)
    hull = Hull(r, z, r0, sigdigits=hull_sigdigits)
    return hull
end

function Hull(hull_data::Dict)
    vs_r = hull_data["vertices_r"]
    vs_z = hull_data["vertices_z"]
    vertices = [[r, z] for (r,z) in zip(vs_r, vs_z)]
    Hull(vertices, hull_data["k"], hull_data["converged"])
end


function Hull(h5_path::String, it_num)
    it_name = @sprintf "iteration_%03d" it_num
    grid_data = h5open(h5_path, "r") do file
        read(file, it_name * "/wind_hull")
    end
    return Hull(grid_data)
end

function Hull(h5_path::String)
    it_keys = h5open(h5_path, "r") do file
        keys(read(file))
    end
    it_nums = [parse(Int, split(key, "_")[end]) for key in it_keys]
    return Hull(h5_path, maximum(it_nums))
end

function is_point_in_wind(hull::ConcaveHull.Hull, point)
    return ConcaveHull.in_hull(point, hull)
end
is_point_in_wind(wi::WindInterpolator, point) = is_point_in_wind(wi.hull, point)
is_point_in_wind(hull::ConcaveHull.Hull, r, z) = is_point_in_wind(hull, [r, z])
