using Distances
export compute_cell_intersection,
    next!, GridIterator, compute_xray_tau, compute_uv_tau, set_iterator!, next_intersection!

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
    rad = b^2 - 4 * a * c
    #println("a $a b $b ri $ri")
    #println("rad $rad")
    if rad < 0
        return NaN
    end
    rad = sqrt(rad)
    sol1 = (-b + rad) / (2a)
    sol2 = (-b - rad) / (2a)
    #println("sol1 $sol1 sol2 $sol2")
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
        #println("index $index")
        #println(length(iterator.r_range))
        if index > length(iterator.r_range)
            #println("BREAKING")
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
        #println("n_steps $n_steps")
        #println("l $lambda")
        #println("r $r")
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
        r = sqrt(x^2+y^2)
        if direction == "r"
            r = r_range[searchsorted_nearest(r_range, r)]
        else
            z = z_range[searchsorted_nearest(z_range, z)]
        end
    end
    return r, z
end

mutable struct GridIterator{T} <: CellIterator{T}
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
    direction_r::Int
    direction_z::Int
    current_r_idx::Int
    current_z_idx::Int
    intersection::Vector{Float64}
    finished::Bool
    function GridIterator(r_range, z_range)
        r_range = round.(r_range, digits = 6)
        z_range = round.(z_range, digits = 6)
        new{typeof(r_range[1])}(
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
    #if (ri * rf * cos(phif-phii) - ri^2 >= 0)
    #    direction_r = 1
    #else
    #    direction_r = -1
    #end
    #(rf == ri) && (direction_r = 1)
    iterator.change_direction_r =
        (ri^2 - ri * rf * cos(phif - phii)) / (ri^2 + rf^2 - 2 * ri * rf * cos(phif - phii))
    update_direction_r!(iterator, 0.0)
    iterator.direction_z = sign(zf - zi)
    current_r_idx = searchsorted_first(iterator.r_range, ri, iterator.direction_r)
    current_z_idx = searchsorted_first(iterator.z_range, zi, iterator.direction_z)
    #println("first $current_r_idx")
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
    #println(iterator.intersection)
    #println(iterator.current_r_idx)
    #println(iterator.r_range)
    #println(iterator.direction_r)
    index =
        min(max(iterator.current_r_idx + iterator.direction_r, 1), length(iterator.r_range))
    return iterator.r_range[index]
end

function next_z(iterator::GridIterator)
    index = iterator.current_z_idx + iterator.direction_z
    return iterator.z_range[index]
end
#next_z(iterator::GridIterator) =
#    iterator.z_range[max(1, abs(iterator.current_z_idx + iterator.direction_z))]
#previous_r(iterator::GridIterator) =
#    iterator.r_range[iterator.current_r_idx - iterator.direction_r]
#previous_z(iterator::GridIterator) =
#    iterator.z_range[iterator.current_z_idx - iterator.direction_z]

function GridIterator(interpolator::WindInterpolator, ri, zi, rf, zf)
    return GridIterator(
        interpolator.density_grid.r_range,
        interpolator.density_grid.z_range,
        ri,
        zi,
        rf,
        zf,
    )
end

function pick_lambda(lambda_r, lambda_z)
    #println("lambda_r $lambda_r")
    #println("lambda_z $lambda_z")
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
    #println("lr $lambda_r")
    #println("lz $lambda_z")
    lambda = pick_lambda(lambda_r, lambda_z)
    iterator.lambda = lambda
    x = iterator.xi * (1 - lambda) + lambda * iterator.xf
    y = iterator.yi * (1 - lambda) + lambda * iterator.yf
    z = iterator.zi * (1 - lambda) + lambda * iterator.zf
    #phi = atan(y, x)
    phi = atan(y, x)
    r = sqrt(x^2 + y^2)
    iterator.intersection[1] = r
    iterator.intersection[2] = phi
    iterator.intersection[3] = z
    #if (iterator.current_r_idx == 1) && (iterator.direction_r == -1)
    #    # time to reverse
    #    iterator.direction_r = 1
    #end

    #epsilon_z = 1e-5 * iterator.direction_z

    #iterator.current_r_idx =
    #    searchsorted_first(iterator.r_range, r+epsilon_r, iterator.direction_r)
    #iterator.current_z_idx =
    #    searchsorted_first(iterator.z_range, z+epsilon_z, iterator.direction_z)

    #epsilon_r = 1e-5 * iterator.direction_r
    #iterator.current_r_idx =
    #    searchsorted_first(iterator.r_range, r+epsilon_r, iterator.direction_r)
    #adjust_r!(iterator, r)
    update_direction_r!(iterator, lambda)
    if lambda_r < lambda_z || isnan(lambda_z)
        iterator.current_r_idx = r_index
    elseif lambda_z < lambda_r || isnan(lambda_r)
        iterator.current_z_idx += iterator.direction_z
    else
        iterator.current_r_idx = r_index
        iterator.current_z_idx += iterator.direction_z
    end
    #println("lr $lambda_r lz $lambda_z ridx $(iterator.current_r_idx)")
end

#function adjust_r!(iterator, r)
#    println("r $r")
#    println("$(iterator.r_range)")
#    rr = current_r(iterator)
#    if iterator.direction_r == 1
#        while r < rr
#            iterator.current_r_idx += 1
#            rr = current_r(iterator)
#        end
#    else
#        while r > rr
#            iterator.current_r_idx -= 1
#            rr = current_r(iterator)
#        end
#    end
#end

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
    #println("\n")
    #println("r0 $r0, r1 $r1, rf $rf")
    #println("z0 $z0, z1 $z1, zf $zf")
    #println(r_cond)
    #println(z_cond)
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
    #println("current r $r0")
    #println("current z $z0")
    #println("finding r intersection to $r1")
    #println("finding z intersection to $z1")
    lambda_r, r_index = get_time_to_intersection_r(iterator)
    lambda_z = get_time_to_intersection_z(z1, iterator.zi, iterator.zf)
    #println("lr $lambda_r lz $lambda_z")
    if is_in_final_cell(iterator, r0, r1, z0, z1, phi0, lambda_r, lambda_z)
        iterator.finished = true
        iterator.intersection[1] = iterator.rf
        iterator.intersection[2] = iterator.phif
        iterator.intersection[3] = iterator.zf
        return
    end
    step_ray!(iterator, lambda_r, lambda_z, r_index)
end

function ionization_cell_xi_kernel(
    xray_luminosity,
    cell_density,
    distance_from_source,
    taux0,
    dx,
    target_ionization_parameter = 1e5,
)
    if dx < distance_from_source
        error()
    else
        ret =
            log(xray_luminosity / (target_ionization_parameter * cell_density * dx^2)) -
            cell_density * SIGMA_T * (dx - distance_from_source) - taux0
        return ret
    end
end

function compute_xray_tau_cell(
    xray_opacity::Boost,
    intersection_size::T,
    distance_from_source::T,
    cell_density::T,
    taux0::T,
    xray_luminosity::T,
    Rg::T;
    atol = 0,
    rtol = 1e-2,
) where {T<:AbstractFloat}
    f(t) = ionization_cell_xi_kernel(
        xray_luminosity,
        cell_density,
        distance_from_source,
        taux0,
        t,
    )
    x1 = distance_from_source
    x2 = distance_from_source + intersection_size
    f1 = f(x1)
    f2 = f(x2)
    if f1 > 0 && f2 > 0
        return taux0 + cell_density * SIGMA_T * intersection_size
    elseif f1 < 0 && f2 < 0
        return taux0 + cell_density * SIGMA_T * intersection_size * 100
    else
        dx0 = find_zero(f, (x1, x2), Bisection())
        return taux0 +
               cell_density *
               SIGMA_T *
               (
                   (dx0 - distance_from_source) +
                   100 * (distance_from_source + intersection_size - dx0)
               )
    end
end

function compute_xray_tau_cell(
    xray_opacity::Thomson,
    intersection_size::T,
    distance_from_source::T,
    cell_density::T,
    taux0::T,
    xray_luminosity::T,
    Rg::T;
    atol = 0,
    rtol = 1e-2,
) where {T<:AbstractFloat}
    return taux0 + intersection_size * cell_density * SIGMA_T
end

function get_density(grid::InterpolationGrid, iterator::GridIterator)
    ridx = iterator.current_r_idx
    r = iterator.r_range[max(iterator.current_r_idx, 1)]
    z = iterator.z_range[iterator.current_z_idx]
    #println("getting density at r $r z $z")
    density = grid.grid[max(iterator.current_r_idx, 1), iterator.current_z_idx]
    #println("density $density")
    return density
end


function dist_to_intersection(iterator, point)
    r0, phi0, z0 = iterator.intersection
    r1, phi1, z1 = point
    # max rquired here because of roundoff errors
    return sqrt(max(0.0, r0^2 + r1^2 + (z1 - z0)^2 - 2.0 * r0 * r1 * cos(phi0 - phi1)))
end

function compute_xray_tau(
    grid::InterpolationGrid{T},
    iterator::GridIterator{T},
    xray_opacity::XRayOpacity,
    ri::T,
    zi::T,
    rf::T,
    zf::T,
    xray_luminosity::T,
    Rg::T,
) where {T<:AbstractFloat}
    if grid.grid === nothing
        return 0.0
    end
    #iterator = grid.iterator
    set_iterator!(iterator, ri, 0.0, zi, rf, 0.0, zf)
    initial_point = [ri, 0.0, zi]
    previous_point = copy(iterator.intersection)
    ret = 0.0
    while !iterator.finished
        previous_point[1] = iterator.intersection[1]
        previous_point[2] = iterator.intersection[2]
        previous_point[3] = iterator.intersection[3]
        cell_density = get_density(grid, iterator)
        #cell_density = get_density(grid, iterator.intersection)
        distance_from_source = dist_to_intersection(iterator, initial_point) * Rg
        next_intersection!(iterator)
        intersection_size = dist_to_intersection(iterator, previous_point) * Rg
        ret = compute_xray_tau_cell(
            xray_opacity,
            intersection_size,
            distance_from_source,
            cell_density,
            ret,
            xray_luminosity,
            Rg,
        )
        #println("previous $(previous_point)")
        #println("next point $(iterator.intersection)")
        #println("---")
        #if ret > 40
        #    return ret
        #end
    end
    return ret
end

function compute_xray_tau(
    grid::InterpolationGrid,
    xray_opacity::XRayOpacity,
    ri,
    zi,
    rf,
    zf,
    xray_luminosity,
    Rg,
)
    compute_xray_tau(
        grid::InterpolationGrid,
        grid.iterator,
        xray_opacity,
        ri,
        zi,
        rf,
        zf,
        xray_luminosity,
        Rg,
    )
end

compute_uv_tau_cell(intersection_size::Float64, cell_density::Float64)::Float64 =
    intersection_size * cell_density

function compute_uv_tau(
    grid::InterpolationGrid,
    iterator::CellIterator,
    ri,
    phii,
    zi,
    rf,
    phif,
    zf,
    Rg,
)
    if grid.grid === nothing
        return 0.0
    end
    #iterator = grid.iterator
    set_iterator!(iterator, ri, phii, zi, rf, phif, zf)
    previous_point = copy(iterator.intersection)
    ret = 0.0
    dist = 0.0
    while !iterator.finished
        previous_point[1] = iterator.intersection[1]
        previous_point[2] = iterator.intersection[2]
        previous_point[3] = iterator.intersection[3]
        cell_density = get_density(grid, iterator)
        next_intersection!(iterator)
        intersection_size = dist_to_intersection(iterator, previous_point)
        #println("previous $(previous_point)")
        #println("next point $(iterator.intersection)")
        #println("index $(iterator.current_r_idx)")
        #println("---")
        dist += intersection_size
        ret += compute_uv_tau_cell(intersection_size, cell_density)
    end
    #println("total distance $dist")
    return ret * Rg * SIGMA_T
end

compute_uv_tau(grid::InterpolationGrid, ri, phii, zi, rf, phif, zf, Rg) =
    compute_uv_tau(grid::InterpolationGrid, grid.iterator, ri, phii, zi, rf, phif, zf, Rg)

compute_uv_tau(grid::InterpolationGrid, ri, phii, rf, zf, Rg) =
    compute_uv_tau(grid::InterpolationGrid, grid.iterator, ri, phii, 0.0, rf, 0.0, zf, Rg)
