using Dates
import Base: flush
export searchsorted_nearest,
    searchsorted_first,
    countsignchanges,
    remove_close_elements,
    iter_paths,
    get_value_in_path,
    set_value_in_path!,
    parse_cl,
    make_cosma_scripts,
    d_euclidean,
    get_time,
    flush

get_time() = Dates.format(now(), "HH:MM:SS")
d_euclidean(r0, r1, z0, z1) = sqrt((r0-r1)^2 + (z0-z1)^2)

function flush()
    flush(stdout)
    flush(stderr)
end

function searchsorted_nearest(a, x)
    idx = searchsortedfirst(a, x)
    if (idx == 1)
        return idx
    end
    if (idx > length(a))
        return length(a)
    end
    if (a[idx] == x)
        return idx
    end
    if (abs(a[idx] - x) < abs(a[idx - 1] - x))
        return idx
    else
        return idx - 1
    end
end

function searchsorted_first(a, x, direction = 1)
    idx = searchsortedfirst(a, x)
    if (idx > length(a))
        return length(a)
    end
    if a[idx] == x
        return idx
    end
    if (idx == 1)
        return 1
    end
    (direction == 0) && (direction = 1)
    return Int(idx - (direction + 1) / 2)
end

function countsignchanges(array::Vector{Float64}, reference = 0)
    counter = 0
    if length(array) == 0
        return 0
    end
    current_sign = sign(array[1] + reference)
    for elem in array
        if sign(elem + reference) != current_sign
            current_sign = sign(elem + reference)
            counter += 1
        end
    end
    return counter
end

#function remove_close_elements(args...; digits = 4)
#    ret = [[array[1] for array in args]]
#    ret = hcat(ret...)
#    for i = 2:length(args[1])
#        noinsert = 0
#        toinsert = zeros(length(args))
#        for (j, array) in enumerate(args)
#            element = trunc(array[i], digits = digits)
#            if element in ret[j, :]
#                noinsert += 1
#            end
#            toinsert[j] = element
#        end
#        if noinsert < length(args)
#            ret = hcat(ret, toinsert)
#        end
#    end
#    return [ret[i, :] for i = 1:length(args)]
#end

function iter_paths(dict)
    function iter(d, path)
        paths = []
        for (k, v) in pairs(d)
            if typeof(v) <: Dict
                paths = vcat(paths, iter(v, vcat(path, k)))
            end
            push!(paths, (vcat(path, k), v))
        end
        return paths
    end
    return iter(dict, [])
end

function get_value_in_path(d, path)
    i = 1
    while i < length(path)
        d = d[path[i]]
        i += 1
    end
    return d[path[end]]
end

function set_value_in_path!(d, path, value)
    i = 1
    while i < length(path)
        if !(path[i] in keys(d))
            d[path[i]] = Dict()
        end
        d = d[path[i]]
        i += 1
    end
    d[path[end]] = value
end

