using LinearAlgebra
using QuadGK

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
    compute_lowest_intercepting_angle(rectangle)

Considering all straight lines passing through the centre (0,0),
returns the lowest theta (measured from the z-axis) of the line
that intersects with the given rectangle.
Note: this assumes that the rectangle is located in the first quadrant.
"""
compute_lowest_intercepting_angle(rectangle) = atan(rectangle.rmin / rectangle.zmax)

"""
    compute_highest_intercepting_angle(rectangle)

Considering all straight lines passing through the centre (0,0),
returns the largest theta (measured from the z-axis) of the line
that intersects with the given rectangle.
Note: this assumes that the rectangle is located in the first quadrant.
"""
compute_highest_intercepting_angle(rectangle) = atan(rectangle.rmax / rectangle.zmin)


"""
    compute_first_intersection(rectangle, theta)

Computes the first point of contact between the rectangle and a line
anchored at the origin with angle (respect to z-axis) theta.
It assumes that they intersect.
"""
function compute_first_intersection(rectangle, theta)
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
    compute_second_intersection(rectangle, theta)

Computes the second point of contact between the rectangle and a line
anchored at the origin with angle (respect to z-axis) theta.
It assumes that they intersect.
"""
function compute_second_intersection(rectangle, theta)
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
    compute_intersection_size(rectangle, theta)

given a line anchored at the origin with angle theta (respect to z-axis),
computes the length of the line segment inside the given rectangle.
"""
function compute_intersection_size(rectangle, theta)
    first_intersection = compute_first_intersection(rectangle, theta)
    second_intersection = compute_second_intersection(rectangle, theta)
    return sqrt(
        (first_intersection[1] - second_intersection[1])^2 +
        (first_intersection[2] - second_intersection[2])^2,
    )
end

"""
    compute_distance_to_first_intersection(rectangle, theta)

given a line anchored at the origin with angle theta (respect to z-axis),
computes the distance from the centre to the first intersection.
"""
function compute_distance_to_first_intersection(rectangle, theta)
    first_intersection = compute_first_intersection(rectangle, theta)
    return norm(first_intersection)
end


"""
    compute_cell_optical_depth(rectangle, theta; density, Rg, mu_electron)

computes the optical depth through a cell
"""
function compute_cell_optical_depth(rectangle, theta; density, Rg, mu_electron = 1.17)
    intersection_size = compute_intersection_size(rectangle, theta)
    return intersection_size * density * SIGMA_T * Rg
end


"""
    compute_optical_depth_to_cell(density_grid; rectangle, theta, Rg, mu_electron)

computes the optical depth to the cell
"""
function compute_optical_depth_to_cell(
    density_grid;
    rectangle,
    theta,
    Rg,
    mu_electron = 1.17,
)
    first_intersection = compute_first_intersection(rectangle, theta)
    return compute_tau_uv(
        density_grid,
        ri = 0.0,
        phii = 0.0,
        zi = 0.0,
        rf = first_intersection[1],
        phif = 0.0,
        zf = first_intersection[2],
        Rg = Rg,
        mu_electron = mu_electron,
        max_tau = Inf,
    )
end


"""
    compute_luminosity_absorbed_by_cell(density_grid; rectangle, Rg, mu_electron, source_luminosity)

Computes the amount of luminosity absorbed by a cell from a central source.
"""
function compute_luminosity_absorbed_by_cell(
    density_grid;
    cell,
    Rg,
    source_luminosity,
    cell_density,
    mu_electron = 1.17,
)
    theta_min = compute_lowest_intercepting_angle(cell)
    theta_max = compute_highest_intercepting_angle(cell)
    function f(theta)
        tau_to_cell = compute_optical_depth_to_cell(
            density_grid,
            rectangle = cell,
            theta = theta,
            Rg = Rg,
            mu_electron = mu_electron,
        )
        tau_cell = compute_cell_optical_depth(
            cell,
            theta,
            density = cell_density,
            Rg = Rg,
            mu_electron = mu_electron,
        )
        return sin(theta) * exp(-tau_to_cell) * (1 - exp(-tau_cell))
    end
    integ, _ = quadgk(f, theta_min, theta_max, atol = 0, rtol = 1e-2)
    return source_luminosity / 2 * integ
end
