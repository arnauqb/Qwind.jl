using LinearAlgebra
using QuadGK

export Rectangle,
    compute_luminosity_absorbed_by_cell,
    compute_luminosity_absorbed_grid,
    compute_scattered_flux_in_cell,
    compute_total_flux_in_cell,
    compute_total_flux_grid

struct Rectangle{T<:AbstractFloat}
    rmin::T
    rmax::T
    zmin::T
    zmax::T
end

function atan_positive(y, x)
    angle = atan(y, x)
    if angle < 0
        return angle + 2π
    else
        return angle
    end
end

function Rectangle(; rmin, rmax, zmin, zmax)
    return Rectangle(Float64(rmin), Float64(rmax), Float64(zmin), Float64(zmax))
end

function change_origin(rectangle, px, py)
    return Rectangle(
        rectangle.rmin - px,
        rectangle.rmax - px,
        rectangle.zmin - py,
        rectangle.zmax - py,
    )
end

"""
    compute_intercepting_angles(rectangle, point)

Considering all straight lines passing through the origin,
returns the lowest and highest angle of the lines that
intercept the rectangle.
"""
function compute_intercepting_angles(rectangle)
    angle1 = atan_positive(rectangle.rmin, rectangle.zmin)
    angle2 = atan_positive(rectangle.rmin, rectangle.zmax)
    angle3 = atan_positive(rectangle.rmax, rectangle.zmin)
    angle4 = atan_positive(rectangle.rmax, rectangle.zmax)
    amin = min(angle1, angle2, angle3, angle4)
    amax = max(angle1, angle2, angle3, angle4)
    return amin, amax
end
function compute_intercepting_angles(rectangle, point)
    new_rect = change_origin(rectangle, point[1], point[2])
    return compute_intercepting_angles(new_rect)
end

"""
    compute_first_intersection(rectangle, theta)

Computes the first point of contact between the rectangle and a line
anchored at the origin with angle (respect to z-axis) theta.
It assumes that they intersect.
"""
function compute_first_intersection(rectangle, theta)
    cotheta = cot(theta)
    # first quadrant
    if 0 <= theta <= π / 2
        zi = cotheta * rectangle.rmin
        if zi < rectangle.zmin
            ri = rectangle.zmin / cotheta
            return [ri, rectangle.zmin]
        else
            return [rectangle.rmin, zi]
        end
        # fourth quadrant
    elseif π / 2 <= theta <= π
        zi = cotheta * rectangle.rmin
        if zi > rectangle.zmax
            ri = rectangle.zmax / cotheta
            return [ri, rectangle.zmax]
        else
            return [rectangle.rmin, zi]
        end
        # third quadrant
    elseif π <= theta <= 3π / 2
        zi = cotheta * rectangle.rmax
        if zi > rectangle.zmax
            ri = rectangle.zmax / cotheta
            return [ri, rectangle.zmax]
        else
            return [rectangle.rmax, zi]
        end
        # second quadrant
    else
        zi = cotheta * rectangle.rmax
        if zi < rectangle.zmin
            ri = rectangle.zmin / cotheta
            return [ri, rectangle.zmin]
        else
            return [rectangle.rmax, zi]
        end
    end
end
function compute_first_intersection(rectangle, theta, point)
    rectangle_centered_point = change_origin(rectangle, point[1], point[end])
    intersection = compute_first_intersection(rectangle_centered_point, theta)
    return intersection[1] + point[1], intersection[2] + point[end]
end

"""
    compute_second_intersection(rectangle, theta)

Computes the second point of contact between the rectangle and a line
anchored at the origin with angle (respect to z-axis) theta.
It assumes that they intersect.
"""
function compute_second_intersection(rectangle, theta)
    cotheta = cot(theta)
    # first quadrant
    if 0 <= theta <= π / 2
        zi = cotheta * rectangle.rmax
        if zi > rectangle.zmax
            ri = rectangle.zmax / cotheta
            return [ri, rectangle.zmax]
        else
            return [rectangle.rmax, zi]
        end
        # fourth quadrant
    elseif π / 2 <= theta <= π
        zi = cotheta * rectangle.rmax
        if zi < rectangle.zmin
            ri = rectangle.zmin / cotheta
            return [ri, rectangle.zmin]
        else
            return [rectangle.rmax, zi]
        end
        # third quadrant
    elseif π <= theta <= 3π / 2
        zi = cotheta * rectangle.rmin
        if zi < rectangle.zmin
            ri = rectangle.zmin / cotheta
            return [ri, rectangle.zmin]
        else
            return [rectangle.rmin, zi]
        end
        # second quadrant
    else
        zi = cotheta * rectangle.rmin
        if zi > rectangle.zmax
            ri = rectangle.zmax / cotheta
            return [ri, rectangle.zmax]
        else
            return [rectangle.rmin, zi]
        end
    end

end
function compute_second_intersection(rectangle, theta, point)
    rectangle_centered_point = change_origin(rectangle, point[1], point[end])
    intersection = compute_second_intersection(rectangle_centered_point, theta)
    return intersection + point[1], intersection + point[end]
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
    origin,
    rectangle,
    theta,
    Rg,
    mu_electron = 1.17,
)
    first_intersection = compute_first_intersection(rectangle, theta, origin)
    return compute_tau_uv(
        density_grid,
        ri = origin[1],
        phii = origin[2],
        zi = origin[3],
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
    source_position,
    cell,
    Rg,
    source_luminosity,
    cell_density,
    mu_electron = 1.17,
)
    cell_centered_source = change_origin(cell, source_position[1], source_position[2])
    theta_min, theta_max = compute_intercepting_angles(cell_centered_source)
    function f(theta)
        tau_to_cell = compute_optical_depth_to_cell(
            density_grid,
            origin = source_position,
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
    integ, _ = quadgk(f, theta_min, theta_max, atol = 0, rtol = 1e-1)
    return source_luminosity / (4π) * integ
end

function compute_luminosity_absorbed_grid(
    dgrid;
    Rg,
    source_luminosity,
    source_position,
    mu_electron = 1.17,
)
    ret = zeros(length(dgrid.r_range) - 1, length(dgrid.z_range) - 1)
    for i = 1:(length(dgrid.r_range) - 1)
        for j = 1:(length(dgrid.z_range) - 1)
            cell = Rectangle(
                rmin = dgrid.r_range[i],
                rmax = dgrid.r_range[i + 1],
                zmin = dgrid.z_range[j],
                zmax = dgrid.z_range[j + 1],
            )
            lumin_absorbed = compute_luminosity_absorbed_by_cell(
                dgrid,
                source_position = source_position,
                cell = cell,
                Rg = Rg,
                cell_density = dgrid.grid[i, j],
                source_luminosity = source_luminosity,
                mu_electron = mu_electron,
            )
            ret[i, j] = lumin_absorbed
        end
    end
    return ret
end



function compute_scattered_flux_in_cell(
    density_grid;
    scattered_luminosity_per_cell,
    cell,
    cell_density,
    Rg,
    mu_electron,
)
    ret = 0.0
    for i = 1:(length(density_grid.r_range) - 1)
        for j = 1:(length(density_grid.z_range) - 1)
            r_source = (density_grid.r_range[i + 1] - density_grid.r_range[i]) / 2
            z_source = (density_grid.z_range[j + 1] - density_grid.z_range[j]) / 2
            source_luminosity = scattered_luminosity_per_cell[i, j]
            cell_density = density_grid.grid[i, j]
            r_target = (cell.rmin + cell.rmax) / 2
            z_target = (cell.zmin + cell.zmax) / 2
            function f(phi)
                distance = compute_distance_cylindrical(
                    r_source,
                    phi,
                    z_source,
                    r_target,
                    0,
                    z_target,
                ) * Rg
                tau = compute_tau_uv(
                    density_grid,
                    ri = r_source,
                    phii = phi,
                    zi = z_source,
                    rf = r_target,
                    phif = 0,
                    zf = z_target,
                    Rg = Rg,
                    mu_electron = mu_electron,
                    max_tau = 50,
                )
                return source_luminosity * exp(-tau) / distance^2
            end
            #f(phi) =
            #f(phi) = compute_luminosity_absorbed_by_cell(
            #    density_grid,
            #    source_position = [r_source, phi, z_source],
            #    source_luminosity = source_luminosity,
            #    cell = cell,
            #    Rg = Rg,
            #    mu_electron = mu_electron,
            #    cell_density = cell_density,
            #)
            flux, _ = 2 .* quadgk(f, 0, π, atol = 0, rtol = 1e-2)
            ret += flux
        end
    end
    return ret
end

function compute_total_flux_in_cell(
    density_grid;
    scattered_luminosity_per_cell,
    source_position,
    source_luminosity,
    cell,
    cell_density,
    Rg,
    mu_electron,
)
    scattered = compute_scattered_flux_in_cell(
        density_grid,
        scattered_luminosity_per_cell = scattered_luminosity_per_cell,
        cell = cell,
        cell_density = cell_density,
        Rg = Rg,
        mu_electron = mu_electron,
    )
    rf = (cell.rmin + cell.rmax) / 2
    zf = (cell.zmin + cell.zmax) / 2
    ridx = searchsortedfirst(density_grid.r_range, cell.rmin)
    zidx = searchsortedfirst(density_grid.z_range, cell.zmin)
    distance_center_sq = (rf^2 + zf^2) * Rg^2
    tau_center = compute_tau_uv(
        density_grid,
        ri = 0,
        phii = 0,
        zi = 0,
        rf = rf,
        phif = 0,
        zf = zf,
        Rg = Rg,
        mu_electron = mu_electron,
        max_tau = 50,
    )
    from_center =
        scattered_luminosity_per_cell[ridx, zidx] * exp(-tau_center) / distance_center_sq
    return from_center + scattered
end

function compute_total_flux_grid(density_grid; Rg, mu_electron = 1.17, source_luminosity)
    ret = zeros(length(density_grid.r_range) - 1, length(density_grid.z_range) - 1)
    absorbed_from_center = compute_luminosity_absorbed_grid(
        density_grid,
        Rg = Rg,
        source_position = [0, 0, 0],
        source_luminosity = source_luminosity,
    )
    function f(i)
        println(i)
        ret = zeros(length(density_grid.z_range) -1)
        for j = 1:(length(density_grid.z_range) - 1)
            cell = Rectangle(
                density_grid.r_range[i],
                density_grid.r_range[i + 1],
                density_grid.z_range[j],
                density_grid.z_range[j + 1],
            )
            ret[j] = compute_total_flux_in_cell(
                density_grid,
                scattered_luminosity_per_cell = absorbed_from_center,
                cell = cell,
                cell_density = density_grid.grid[i, j],
                Rg = Rg,
                mu_electron = mu_electron,
                source_position = [0, 0, 0],
                source_luminosity = source_luminosity,
            )
        end
        return ret
    end
    rets = pmap(f, 1:length(density_grid.r_range)-1)
    rets = reduce(hcat, rets)'
    return rets
end
