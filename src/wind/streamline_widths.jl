export KDTree, get_streamlines_widths

using NearestNeighbors

function NearestNeighbors.KDTree(sl::Streamline)
    points = hcat(sl.r, sl.z)'
    return KDTree(points)
end
function NearestNeighbors.KDTree(sl::Sundials.IDAIntegrator)
    points = hcat(sl.p.data[:r], sl.p.data[:z])'
    return KDTree(points)
end

function get_nearest_point(kdtree, point)
    idxs, dists = knn(kdtree, point, 1)
    return kdtree.data[idxs[1]]
end

function get_distance_to_line(kdtree, point)
    idcs, dists = knn(kdtree, point, 1)
    return dists[1]
end

function get_closest_point_and_distance(kdtree, point)
    idcs, dists = knn(kdtree, point, 1)
    return kdtree.data[idcs[1]], dists[1]
end

function get_distances_between_lines(line1, line2)
    ret = zero(line1.r)
    kdt2 = KDTree(line2)
    for i in 1:length(line1.r)
        r = line1.r[i]
        z = line1.z[i]
        point = [r, z]
        vr = line1.vr[i]
        vz = line1.vz[i]
        v = sqrt(vr^2 + vz^2)
        closest_point, distance = get_closest_point_and_distance(kdt2, point)
        p_vector = closest_point - point
        cosθ = (p_vector[1] * vr + p_vector[2] * vz) / (v * distance)
        #ret[i] = abs(distance / cosθ)
        ret[i] = distance
    end
    return ret
end

function get_distances_between_lines(
    line1::Sundials.IDAIntegrator,
    line2::Sundials.IDAIntegrator,
)
    ret = zero(line1.p.data[:r])
    kdt2 = KDTree(line2)
    for i in 1:length(line1.p.data[:r])
        r = line1.p.data[:r][i]
        z = line1.p.data[:z][i]
        point = [r, z]
        vr = line1.p.data[:vr][i]
        vz = line1.p.data[:vz][i]
        v = sqrt(vr^2 + vz^2)
        closest_point, distance = get_closest_point_and_distance(kdt2, point)
        p_vector = closest_point - point
        cosθ = (p_vector[1] * vr + p_vector[2] * vz) / (v * distance)
        ret[i] = abs(distance / cosθ)
    end
    return ret
end

function get_streamlines_widths(streamlines)
    ret = []
    for i = 1:(length(streamlines) - 1)
        push!(ret, get_distances_between_lines(streamlines[i], streamlines[i + 1]))
    end
    return ret
end
