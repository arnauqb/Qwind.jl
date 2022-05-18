export KDTree, get_line_widths

using NearestNeighbors

function NearestNeighbors.KDTree(sl::Streamline)
    points = hcat(sl.r, sl.z)'
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

function get_line_widths(line::Streamline, kdt1, kdt2)
    ret = zero(line.r)
    for (i, (r, z)) in enumerate(zip(line.r, line.z))
        point = [r, z]
        ret[i] = get_distance_to_line(kdt1, point) + get_distance_to_line(kdt2, point)
    end
    return ret
end
