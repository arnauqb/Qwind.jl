using LinearAlgebra

export Rectangle

struct Rectangle{T<:AbstractFloat}
    rmin::T
    rmax::T
    zmin::T
    zmax::T
end

function Rectangle(; rmin, rmax, zmin, zmax)
    return Rectangle(Float64(rmin), Float64(rmax), Float64(zmin), Float64(zmax))
end


"""
    get_lowest_intercepting_angle(rectangle)

Considering all straight lines passing through the centre (0,0),
returns the lowest theta (measured from the z-axis) of the line
that intersects with the given rectangle.
Note: this assumes that the rectangle is located in the first quadrant.
"""
get_lowest_intercepting_angle(rectangle) = atan(rectangle.rmin / rectangle.zmax)

"""
    get_highest_intercepting_angle(rectangle)

Considering all straight lines passing through the centre (0,0),
returns the largest theta (measured from the z-axis) of the line
that intersects with the given rectangle.
Note: this assumes that the rectangle is located in the first quadrant.
"""
get_highest_intercepting_angle(rectangle) = atan(rectangle.rmax / rectangle.zmin)


"""
    get_first_intersection(rectangle, theta)

Computes the first point of contact between the rectangle and a line
anchored at the origin with angle (respect to z-axis) theta.
It assumes that they intersect.
"""
function get_first_intersection(rectangle, theta)
    cotheta = cot(theta)
    zi = cotheta * rectangle.rmin
    if zi < rectangle.zmin
        ri = rectangle.zmin / cotheta
        return [ri, rectangle.zmin]
    else
        return [rectangle.rmin, zi]
    end
end

"""
    get_second_intersection(rectangle, theta)

Computes the second point of contact between the rectangle and a line
anchored at the origin with angle (respect to z-axis) theta.
It assumes that they intersect.
"""
function get_second_intersection(rectangle, theta)
    cotheta = cot(theta)
    zi = cotheta * rectangle.rmax
    if zi > rectangle.zmax
        ri = rectangle.zmax / cotheta
        return [ri, rectangle.zmax]
    else
        return [rectangle.rmax, zi]
    end
end

"""
    get_intersection_size(rectangle, theta)

given a line anchored at the origin with angle theta (respect to z-axis),
computes the length of the line segment inside the given rectangle.
"""
function get_intersection_size(rectangle, theta)
    first_intersection = get_first_intersection(rectangle, theta)
    second_intersection = get_second_intersection(rectangle, theta)
    return sqrt(
        (first_intersection[1] - second_intersection[1])^2 +
        (first_intersection[2] - second_intersection[2])^2,
    )
end

"""
    get_distance_to_first_intersection(rectangle, theta)

given a line anchored at the origin with angle theta (respect to z-axis),
computes the distance from the centre to the first intersection.
"""
function get_distance_to_first_intersection(rectangle, theta)
    first_intersection = get_first_intersection(rectangle, theta)
    return norm(first_intersection)
end
