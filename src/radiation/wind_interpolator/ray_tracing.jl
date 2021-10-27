using Distances
export compute_cell_intersection,
    next!, GridIterator, compute_tau_xray, compute_tau_uv, set_iterator!, next_intersection!

function point_outside_grid(r_range, z_range, r, z)
    if !(r_range[1] <= r <= r_range[end]) || !(z_range[1] <= z <= z_range[end])
        return true
    else
        return false
    end
end

function get_time_to_intersection_r(r, a, b, ri, current_lambda)
    if a == 0.0
        return Inf
    end
    c = ri^2 - r^2
    rad = b^2 - 4.0 * a * c
    if rad < 0.0
        return NaN
    end
    rad = sqrt(rad)
    sol1 = (-b + rad) / (2a)
    sol2 = (-b - rad) / (2a)
    if current_lambda > max(sol1, sol2)
        return Inf
    end
    if sol1 > 0
        if sol2 > 0
            if sol1 < sol2
                if sol1 > current_lambda
                    return sol1
                else
                    return sol2
                end
            else
                if sol2 > current_lambda
                    return sol2
                else
                    return sol1
                end
            end
        else
            return sol1
        end
    else
        return sol2
    end
end

function get_time_to_intersection_r(iterator)
    lambda = Inf
    index = 0
    n_steps = 0
    while true
        index = max(abs(iterator.current_r_idx + (n_steps + 1) * iterator.direction_r), 1)
        if index > length(iterator.r_range)
            return Inf, 0.0
        end
        r = iterator.r_range[index]
        lambda = get_time_to_intersection_r(
            r,
            iterator.a,
            iterator.b,
            iterator.ri,
            iterator.lambda,
        )
        n_steps += 1
        if !isnan(lambda)
            break
        end
    end
    return lambda, index
end


function get_time_to_intersection_z(z, zi, zf)
    abs((z - zi) / (zf - zi))
end

function get_intersection_with_grid(
    r,
    z,
    r_range,
    z_range,
    ri,
    phii,
    zi,
    rf,
    phif,
    zf;
    from_outside,
)
    if from_outside
        lambda_r = -Inf
        lambda_z = -Inf
    else
        lambda_r = Inf
        lambda_z = Inf
    end
    if r < r_range[1]
        a = ri^2 + rf^2 - 2 * ri * rf * cos(phif - phii)
        b = 2 * ri * rf * cos(phif - phii) - 2 * ri^2
        lambda_r = get_time_to_intersection_r(r_range[1], a, b, ri, 0.0)
    elseif r > r_range[end]
        a = ri^2 + rf^2 - 2 * ri * rf * cos(phif - phii)
        b = 2 * ri * rf * cos(phif - phii) - 2 * ri^2
        lambda_r = get_time_to_intersection_r(r_range[end], a, b, ri, 0.0)
    end
    if z < z_range[1]
        lambda_z = get_time_to_intersection_z(z_range[1], zi, zf)
    elseif z > z_range[end]
        lambda_z = get_time_to_intersection_z(z_range[end], zi, zf)
    end
    if from_outside
        if lambda_r > lambda_z
            direction = "r"
            lambda = lambda_r
        else
            direction = "z"
            lambda = lambda_z
        end
        #lambda = max(lambda_r, lambda_z)
    else
        #lambda = min(lambda_r, lambda_z)
        if lambda_r < lambda_z
            direction = "r"
            lambda = lambda_r
        else
            direction = "z"
            lambda = lambda_z
        end
    end
    if abs(lambda) != Inf && !isnan(lambda)
        xi = ri * cos(phii)
        yi = ri * sin(phii)
        xf = rf * cos(phif)
        yf = rf * sin(phif)
        x = xi * (1 - lambda) + lambda * xf
        y = yi * (1 - lambda) + lambda * yf
        z = zi * (1 - lambda) + lambda * zf
        r = sqrt(x^2 + y^2)
        if direction == "r"
            r = r_range[searchsorted_nearest(r_range, r)]
        else
            z = z_range[searchsorted_nearest(z_range, z)]
        end
    end
    return r, z
end

mutable struct GridIterator{T,U,V}
    r_range::Vector{T}
    z_range::Vector{T}
    ri::T
    phii::T
    zi::T
    rf::T
    phif::T
    zf::T
    xi::T
    xf::T
    yi::T
    yf::T
    a::T
    b::T
    lambda::T
    change_direction_r::T
    direction_r::U
    direction_z::U
    current_r_idx::U
    current_z_idx::U
    intersection::Vector{T}
    previous_point::Vector{T}
    finished::V
    function GridIterator(r_range, z_range)
        r_range = round.(r_range, digits = 6)
        z_range = round.(z_range, digits = 6)
        new{typeof(r_range[1]),Int,Bool}(
            r_range,
            z_range,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            1,
            1,
            1,
            1,
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            false,
        )
    end
end

function GridIterator(r_range, z_range, ri, phii, zi, rf, phif, zf)
    gi = GridIterator(r_range, z_range)
    set_iterator!(gi, ri, phii, zi, rf, phif, zf)
    return gi
end

function GridIterator(r_range, z_range, ri, zi, rf, zf)
    gi = GridIterator(r_range, z_range)
    set_iterator!(gi, ri, 0.0, zi, rf, 0.0, zf)
    return gi
end


function corner_sign(xcorner, ycorner, x1, y1, x2, y2)
    return ((y2 - y1) * xcorner + (x1 - x2) * ycorner + (x2 * y1 - x1 * y2)) >= 0
end

function intersects_rectangle(r_min, r_max, z_min, z_max, ri, zi, rf, zf)
    # if first or last point inside rectangle then true
    if r_min <= ri <= r_max && z_min <= zi <= z_max
        return true
    elseif r_min <= rf <= r_max && z_min <= zf <= z_max
        return true
    elseif (ri <= r_min && rf <= r_min) || (ri >= r_max && rf >= r_max)
        return false
    elseif (zi <= z_min && zf <= z_min) || (zi >= z_max && rf >= z_max)
        return false
    else
        # here we know that both points are outside the grid, we need to work out
        # if they intersect the rectangle. One way to do it is to check if the segment
        # "lies" in the same side for all corners
        sign1 = corner_sign(r_min, z_min, ri, zi, rf, zf)
        sign2 = corner_sign(r_min, z_max, ri, zi, rf, zf)
        (sign2 != sign1) && return true
        sign3 = corner_sign(r_max, z_min, ri, zi, rf, zf)
        (sign3 != sign1) && return true
        sign4 = corner_sign(r_max, z_max, ri, zi, rf, zf)
        (sign4 != sign1) && return true
    end
    return false
end

function set_iterator!(iterator::GridIterator, ri, phii, zi, rf, phif, zf)
    if !intersects_rectangle(
        iterator.r_range[1],
        iterator.r_range[end],
        iterator.z_range[1],
        iterator.z_range[end],
        ri,
        zi,
        rf,
        zf,
    )
        iterator.ri = ri
        iterator.phii = phii
        iterator.zi = zi
        iterator.rf = rf
        iterator.phif = phif
        iterator.zf = zf
        iterator.finished = true
        iterator.lambda = 0.0
        return
    end
    ri, zi = get_intersection_with_grid(
        ri,
        zi,
        iterator.r_range,
        iterator.z_range,
        ri,
        phii,
        zi,
        rf,
        phif,
        zf,
        from_outside = true,
    )
    rf, zf = get_intersection_with_grid(
        rf,
        zf,
        iterator.r_range,
        iterator.z_range,
        ri,
        phii,
        zi,
        rf,
        phif,
        zf,
        from_outside = false,
    )
    iterator.change_direction_r =
        (ri^2 - ri * rf * cos(phif - phii)) / (ri^2 + rf^2 - 2 * ri * rf * cos(phif - phii))
    update_direction_r!(iterator, 0.0)
    iterator.direction_z = sign(zf - zi)
    current_r_idx = searchsorted_first(iterator.r_range, ri, iterator.direction_r)
    current_z_idx = searchsorted_first(iterator.z_range, zi, iterator.direction_z)
    iterator.ri = ri
    iterator.phii = phii
    iterator.zi = zi
    iterator.rf = rf
    iterator.phif = phif
    iterator.zf = zf
    iterator.a = ri^2 + rf^2 - 2 * ri * rf * cos(phif - phii)
    iterator.b = 2 * ri * rf * cos(phif - phii) - 2 * ri^2
    iterator.xi = iterator.ri * cos(iterator.phii)
    iterator.xf = iterator.rf * cos(iterator.phif)
    iterator.yi = iterator.ri * sin(iterator.phii)
    iterator.yf = iterator.rf * sin(iterator.phif)
    iterator.lambda = 0.0
    iterator.current_r_idx = current_r_idx
    iterator.current_z_idx = current_z_idx
    iterator.intersection[1] = ri
    iterator.intersection[2] = phii
    iterator.intersection[3] = zi
    iterator.previous_point .= iterator.intersection
    iterator.finished = false
    return
end


current_r(iterator::GridIterator) = iterator.r_range[iterator.current_r_idx]
current_z(iterator::GridIterator) = iterator.z_range[iterator.current_z_idx]

function update_direction_r!(iterator, lambda)
    if lambda < iterator.change_direction_r
        iterator.direction_r = -1
    else
        iterator.direction_r = 1
    end
end

function next_r(iterator::GridIterator)
    index =
        min(max(iterator.current_r_idx + iterator.direction_r, 1), length(iterator.r_range))
    return iterator.r_range[index]
end

function next_z(iterator::GridIterator)
    index = iterator.current_z_idx + iterator.direction_z
    return iterator.z_range[index]
end

function pick_lambda(lambda_r, lambda_z)
    if lambda_r < lambda_z || isnan(lambda_z)
        lambda = lambda_r
    elseif lambda_z < lambda_r || isnan(lambda_r)
        lambda = lambda_z
    else
        lambda = lambda_r
    end
    return lambda
end

function step_ray!(iterator::GridIterator, lambda_r, lambda_z, r_index)
    lambda = pick_lambda(lambda_r, lambda_z)
    iterator.lambda = lambda
    x = iterator.xi * (1 - lambda) + lambda * iterator.xf
    y = iterator.yi * (1 - lambda) + lambda * iterator.yf
    z = iterator.zi * (1 - lambda) + lambda * iterator.zf
    phi = atan(y, x)
    r = sqrt(x^2 + y^2)
    iterator.intersection[1] = r
    iterator.intersection[2] = phi
    iterator.intersection[3] = z
    update_direction_r!(iterator, lambda)
    if lambda_r < lambda_z || isnan(lambda_z)
        iterator.current_r_idx = r_index
    elseif lambda_z < lambda_r || isnan(lambda_r)
        iterator.current_z_idx += iterator.direction_z
    else
        iterator.current_r_idx = r_index
        iterator.current_z_idx += iterator.direction_z
    end
end

function is_in_final_cell(iterator, r0, r1, z0, z1, phi0, lambda_r, lambda_z)
    r_cond = false
    z_cond = false
    rf = iterator.rf
    phif = iterator.phif
    zf = iterator.zf
    if isnan(lambda_r) || lambda_r > 1
        r_cond = true
    else
        if r1 >= r0
            if r0 <= rf <= r1
                r_cond = true
            end
        else
            if r0 >= rf >= r1
                r_cond = true
            end
        end
    end
    if isnan(lambda_z) || lambda_z > 1
        z_cond = true
    else
        if z1 >= z0
            if z0 <= zf <= z1
                z_cond = true
            end
        else
            if z0 >= zf >= z1
                z_cond = true
            end
        end
    end
    return r_cond && z_cond
end


function next_intersection!(iterator::GridIterator)
    if iterator.finished
        return nothing
    end
    epsilon = 1e-10
    r0 = iterator.intersection[1]
    phi0 = iterator.intersection[2]
    z0 = iterator.intersection[3]
    r1 = next_r(iterator)
    z1 = next_z(iterator)
    lambda_r, r_index = get_time_to_intersection_r(iterator)
    lambda_z = get_time_to_intersection_z(z1, iterator.zi, iterator.zf)
    if is_in_final_cell(iterator, r0, r1, z0, z1, phi0, lambda_r, lambda_z)
        iterator.finished = true
        iterator.intersection[1] = iterator.rf
        iterator.intersection[2] = iterator.phif
        iterator.intersection[3] = iterator.zf
        return
    end
    step_ray!(iterator, lambda_r, lambda_z, r_index)
end

function ionization_cell_xi_kernel(;
    xray_luminosity,
    cell_density,
    distance_from_source,
    taux0,
    dx,
    mu_nucleon,
    mu_electron,
    target_ionization_parameter = 1e5,
)
    if dx < distance_from_source
        error()
    else
        ret =
            log(
                xray_luminosity /
                (target_ionization_parameter * cell_density / mu_nucleon * dx^2),
            ) - cell_density * mu_electron * SIGMA_T * (dx - distance_from_source) - taux0
        return ret
    end
end

function compute_tau_xray_cell(
    xray_opacity::Boost,
    intersection_size,
    distance_from_source,
    cell_density,
    taux0,
    xray_luminosity,
    Rg,
    mu_nucleon,
    mu_electron;
    atol = 0,
    rtol = 1e-2,
) where {T<:AbstractFloat}
    f(t) = ionization_cell_xi_kernel(
        xray_luminosity = xray_luminosity,
        cell_density = cell_density,
        distance_from_source = distance_from_source,
        taux0 = taux0,
        mu_nucleon = mu_nucleon,
        mu_electron = mu_electron,
        dx = t,
    )
    x1 = distance_from_source
    x2 = distance_from_source + intersection_size
    f1 = f(x1)
    f2 = f(x2)
    if f1 > 0 && f2 > 0
        return taux0 + cell_density * SIGMA_T * mu_electron * intersection_size
    elseif f1 < 0 && f2 < 0
        return taux0 + cell_density * SIGMA_T * mu_electron * intersection_size * 100
    else
        dx0 = find_zero(f, (x1, x2), Bisection())
        return taux0 +
               cell_density *
               mu_electron *
               SIGMA_T *
               (
                   (dx0 - distance_from_source) +
                   100 * (distance_from_source + intersection_size - dx0)
               )
    end
end

function compute_tau_xray_cell(
    xray_opacity::Thomson,
    intersection_size,
    distance_from_source,
    cell_density,
    taux0,
    xray_luminosity,
    Rg,
    mu_nucleon,
    mu_electron;
    atol = 0,
    rtol = 1e-2,
) where {T<:AbstractFloat}
    return taux0 + intersection_size * mu_electron * cell_density * SIGMA_T
end

function get_density(grid::InterpolationGrid, iterator::GridIterator)
    ridx = iterator.current_r_idx
    r = iterator.r_range[max(iterator.current_r_idx, 1)]
    z = iterator.z_range[iterator.current_z_idx]
    density = grid.grid[max(iterator.current_r_idx, 1), iterator.current_z_idx]
    return density
end


function dist_to_intersection(iterator, point)
    r0, phi0, z0 = iterator.intersection
    r1, phi1, z1 = point
    # max rquired here because of roundoff errors
    return sqrt(max(0.0, r0^2 + r1^2 + (z1 - z0)^2 - 2.0 * r0 * r1 * cos(phi0 - phi1)))
end

function compute_tau_xray(
    grid::InterpolationGrid,
    iterator::GridIterator,
    xray_opacity::XRayOpacityFlag,
    ri,
    zi,
    rf,
    zf,
    xray_luminosity,
    Rg,
    mu_nucleon,
    mu_electron,
)
    if grid.grid === nothing
        return 0.0
    end
    set_iterator!(iterator, ri, 0.0, zi, rf, 0.0, zf)
    initial_point = [ri, 0.0, zi]
    ret = 0.0
    while !iterator.finished
        iterator.previous_point .= iterator.intersection
        cell_density = get_density(grid, iterator)
        distance_from_source = dist_to_intersection(iterator, initial_point) * Rg
        next_intersection!(iterator)
        intersection_size = dist_to_intersection(iterator, iterator.previous_point) * Rg
        ret = compute_tau_xray_cell(
            xray_opacity,
            intersection_size,
            distance_from_source,
            cell_density,
            ret,
            xray_luminosity,
            Rg,
            mu_nucleon,
            mu_electron,
        )
    end
    return ret
end

function compute_tau_xray(
    grid::InterpolationGrid,
    xray_opacity::XRayOpacityFlag;
    ri,
    zi,
    rf,
    zf,
    xray_luminosity,
    Rg,
    mu_nucleon,
    mu_electron,
)
    compute_tau_xray(
        grid::InterpolationGrid,
        grid.iterator,
        xray_opacity,
        ri,
        zi,
        rf,
        zf,
        xray_luminosity,
        Rg,
        mu_nucleon,
        mu_electron,
    )
end

compute_tau_uv_cell(intersection_size, cell_density, mu_electron) =
    intersection_size * cell_density * mu_electron

function compute_tau_uv(
    grid::InterpolationGrid,
    iterator::GridIterator,
    ri,
    phii,
    zi,
    rf,
    phif,
    zf,
    Rg,
    mu_electron;
    max_tau = 50.0,
)
    if grid.grid === nothing
        return 0.0
    end
    set_iterator!(iterator, ri, phii, zi, rf, phif, zf)
    ret = 0.0
    sigmarg = Rg * SIGMA_T
    while !iterator.finished
        iterator.previous_point .= iterator.intersection
        cell_density = get_density(grid, iterator)
        next_intersection!(iterator)
        intersection_size = dist_to_intersection(iterator, iterator.previous_point)
        ret += compute_tau_uv_cell(intersection_size, cell_density, mu_electron)
        if ret * sigmarg > max_tau
            return max_tau
        end
    end
    return ret * sigmarg
end

compute_tau_uv(
    grid::InterpolationGrid;
    ri,
    phii,
    zi,
    rf,
    phif,
    zf,
    Rg,
    mu_electron,
    max_tau = 50.0,
) = compute_tau_uv(
    grid::InterpolationGrid,
    grid.iterator,
    ri,
    phii,
    zi,
    rf,
    phif,
    zf,
    Rg,
    mu_electron,
    max_tau = max_tau,
)

