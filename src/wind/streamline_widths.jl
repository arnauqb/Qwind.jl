export KDTree, get_streamlines_widths

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

function get_distances_between_lines(line1, line2)
    ret = zero(line1.r)
    kdt2 = KDTree(line2)
    for (i, (r, z)) in enumerate(zip(line1.r, line1.z))
        point = [r, z]
        ret[i] = get_distance_to_line(kdt2, point)
    end
    return ret
end

function get_streamlines_widths(streamlines)
    ret = []
    for i in 1:length(streamlines)-1
        push!(ret, get_distances_between_lines(streamlines[i], streamlines[i+1]))
    end
    return ret
end
