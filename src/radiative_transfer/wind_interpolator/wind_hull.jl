import ConcaveHull
import ConcaveHull: Hull
export DenseIntegrator, DenseIntegrators, Hull

function filter_array(a, b, value)
    mask = a .!= value
    return a[mask], b[mask]
end

function Hull(r::Vector{Float64}, z::Vector{Float64}, r0::Vector{Float64}; sigdigits=6)
    # remove points that are too close to each other
    points = hcat(r, z)
    points = round.(points, sigdigits = sigdigits)
    points = unique(points, dims = 1)
    points = [[points[i, 1], points[i, 2]] for i = 1:size(points)[1]]
    hull = ConcaveHull.concave_hull(points, 1)
    return hull#, points
end

function Hull(integrators::Vector{<:Sundials.IDAIntegrator}, max_times)
    r0 = [integ.p.r0 for integ in integrators]
    integrators_interpolated_linear = interpolate_integrators(
        integrators,
        max_times = max_times,
        n_timesteps = 100,
        log = true,
    )
    r, z, _, _, _ = reduce_integrators(integrators_interpolated_linear)
    @info "Constructing wind hull"
    hull = nothing
    for sigd in [6, 5, 4]
        hull = Hull(r, z, r0, sigdigits=sigd)
        if !hull.converged
            @warn "Hull did not converge, trying with less sigdigits..."
        end
    end
    if !hull.converged
        error("Hull did not converge!")
    end
    @info "Done"
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
