import ConcaveHull
export construct_wind_hull

function construct_wind_hull(r::Vector{Float64}, z::Vector{Float64}, r0::Vector{Float64}; sigdigits=6)
    # this stores the maximum heights of each streamline, roughly
    z_max = zeros(length(r0))
    for (rp, zp) in zip(r, z)
        ridx = searchsorted_nearest(r0, rp)
        z_max[ridx] = max(z_max[ridx], zp)
    end
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
    return hull
end

function construct_wind_hull(integrators; hull_sigdigits=6)
    r0 = [integ.p.r0 for integ in integrators]
    hull = nothing
    r, z, _, _, _ = reduce_integrators(integrators, no_interpolation=true)
    hull = construct_wind_hull(r, z, r0, sigdigits=hull_sigdigits)
    return hull
end

function is_point_in_wind(hull::ConcaveHull.Hull, point)
    return ConcaveHull.in_hull(point, hull)
end
is_point_in_wind(wi::WindInterpolator, point) = is_point_in_wind(wi.hull, point)
is_point_in_wind(hull::ConcaveHull.Hull, r, z) = is_point_in_wind(hull, [r, z])
