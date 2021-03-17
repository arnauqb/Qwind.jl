import ConcaveHull
export construct_wind_hull

function construct_wind_hull(r::Vector{Float64}, z::Vector{Float64}, r0::Vector{Float64}; sigdigits=2)
    #r = r[z .> 0]
    #z = z[z .> 0]
    #r_log = log10.(r)
    #z_log = log10.(z)
    #r0_log = log10.(r0)
    #r_range_log = log10.(range(minimum(r), maximum(r), step = 1))
    #r_range_log = sort(vcat(r_range_log, r0_log))
    ##z_range_log = range(max(minimum(z_log), -8), maximum(z_log), length = 250)
    #z_range_log = log10.(range(1e-6, maximum(z), step=0.5))

    #zmax_hull = -Inf .* ones(length(r_range_log))
    #rmax_hull = -Inf .* ones(length(z_range_log))
    #rmin_hull = Inf .* ones(length(z_range_log))

    ## top part of the wind
    #for (rp_log, zp_log) in zip(r_log, z_log)
    #    ridx = searchsorted_nearest(r_range_log, rp_log)
    #    zidx = searchsorted_nearest(z_range_log, zp_log)

    #    # get highest position of the wind
    #    zmax_hull[ridx] = max(zmax_hull[ridx], zp_log)

    #    # get furthest out part of the wind
    #    rmax_hull[zidx] = max(rmax_hull[zidx], rp_log)

    #    # get innermost out part of the wind
    #    rmin_hull[zidx] = min(rmin_hull[zidx], rp_log)
    #end

    #r_hull = r0_log
    #z_hull = -10 .* ones(length(r_hull))

    ## filter unassigned
    #mask = zmax_hull .!= -Inf
    #r_hull = vcat(r_hull, r_range_log[mask])
    #z_hull = vcat(z_hull, zmax_hull[mask])

    #mask = rmax_hull .!= -Inf
    #r_hull = vcat(r_hull, rmax_hull[mask])
    #z_hull = vcat(z_hull, z_range_log[mask])

    #mask = rmin_hull .!= Inf
    #r_hull = vcat(r_hull, rmin_hull[mask])
    #z_hull = vcat(z_hull, z_range_log[mask])

    #points = []
    #for (rp, zp) in zip(r_hull, z_hull)
    #    push!(points, [rp, zp])
    #end
    #points = reduce(hcat, points)
    #points = round.(points, sigdigits = sigdigits)
    #points = unique(points, dims = 2)
    #points = [[points[1, i], points[2, i]] for i = 1:size(points)[2]]


    ## new ##
    #points = [[r,z] for (r,z) in zip(r, z)]
    #points = hcat(r_log, z_log)
    ##points = reduce(hcat, points)
    #points = round.(points, sigdigits = sigdigits)
    #points = unique(points, dims = 1)
    #points = [[points[i, 1], points[i, 2]] for i = 1:size(points)[1]]
    #return points
    #
    #
    #
    #new new
    r_range = range(minimum(r), maximum(r), step = 0.1)
    r_range = sort(vcat(r_range, r0))
    z_range = range(0, maximum(z), step=1)

    zmax_hull = -Inf .* ones(length(r_range))
    rmax_hull = -Inf .* ones(length(z_range))
    rmin_hull = Inf .* ones(length(z_range))

    ## top part of the wind
    for (rp, zp) in zip(r, z)
        ridx = searchsorted_nearest(r_range, rp)
        zidx = searchsorted_nearest(z_range, zp)

        # get highest position of the wind
        zmax_hull[ridx] = max(zmax_hull[ridx], zp)

        # get furthest out part of the wind
        rmax_hull[zidx] = max(rmax_hull[zidx], rp)

        # get innermost out part of the wind
        rmin_hull[zidx] = min(rmin_hull[zidx], rp)
    end

    r_hull = r0
    z_hull = zeros(length(r_hull))

    ## filter unassigned
    mask = zmax_hull .!= -Inf
    r_hull = vcat(r_hull, r_range[mask])
    z_hull = vcat(z_hull, zmax_hull[mask])

    mask = rmax_hull .!= -Inf
    r_hull = vcat(r_hull, rmax_hull[mask])
    z_hull = vcat(z_hull, z_range[mask])

    mask = rmin_hull .!= Inf
    r_hull = vcat(r_hull, rmin_hull[mask])
    z_hull = vcat(z_hull, z_range[mask])

    points = [[rp, zp] for (rp, zp) in zip(r,z)]

    #points = [[rp, zp] for (rp, zp) in zip(r_hull, z_hull)]

    points = reduce(hcat, points)
    points = round.(points, sigdigits = sigdigits)
    points = unique(points, dims = 2)
    points = [[points[1, i], points[2, i]] for i = 1:size(points)[2]]
    #####
    
    hull = ConcaveHull.concave_hull(points)
    return hull, points
end

function construct_wind_hull(integrators)
    r0 = [integ.p.r0 for integ in integrators]
    hull = nothing
    r, z, _, _, _ = reduce_integrators(integrators, n_timesteps=500, log=true)
    #dense_integrators = DenseIntegrators(integrators, n_timesteps=100, log=true)
    #r = dense_integrators.r
    #z = dense_integrators.z
    for sigdigits in [6, 5, 4]
        @info "Trying wind hull with $sigdigits sig digits..."
        flush()
        hull = construct_wind_hull(r, z, r0, sigdigits=sigdigits)
        hull.converged && break
    end
    if hull === nothing 
        error("Cannot construct wind hull!")
    end
    return hull
end

function is_point_in_wind(hull::ConcaveHull.Hull, point)
    #return ConcaveHull.in_hull(log10.(point), hull)
    return ConcaveHull.in_hull(point, hull)
end
is_point_in_wind(wi::WindInterpolator, point) = is_point_in_wind(wi.hull, point)
is_point_in_wind(hull::ConcaveHull.Hull, r, z) = is_point_in_wind(hull, [r, z])
