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

get_time_to_intersection(x, xi, xf) = abs((x - xi) / (xf - xi))

function get_intersection_with_grid(r, z, r_range, z_range, ri, rf, zi, zf; from_outside)
    if from_outside
        lambda_r = -Inf
        lambda_z = -Inf
    else
        lambda_r = Inf
        lambda_z = Inf
    end
    if r < r_range[1]
        lambda_r = get_time_to_intersection(r_range[1], ri, rf)
    elseif r > r_range[end]
        lambda_r = get_time_to_intersection(r_range[end], ri, rf)
    end
    if z < z_range[1]
        lambda_z = get_time_to_intersection(z_range[1], zi, zf)
    elseif z > z_range[end]
        lambda_z = get_time_to_intersection(z_range[end], zi, zf)
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
        r = ri + lambda * (rf - ri)
        z = zi + lambda * (zf - zi)
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
    zi::T
    rf::T
    zf::T
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
            1,
            1,
            1,
            1,
            [0.0, 0.0],
            false,
        )
    end
end

function GridIterator(r_range, z_range, ri, zi, rf, zf)
    gi = GridIterator(r_range, z_range)
    set_iterator!(gi, ri, zi, rf, zf)
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
        #Println("sign1 $sign1")
        sign2 = corner_sign(r_min, z_max, ri, zi, rf, zf)
        #println("sign2 $sign2")
        (sign2 != sign1) && return true
        sign3 = corner_sign(r_max, z_min, ri, zi, rf, zf)
        #println("sign3 $sign3")
        (sign3 != sign1) && return true
        sign4 = corner_sign(r_max, z_max, ri, zi, rf, zf)
        #println("sign4 $sign4")
        (sign4 != sign1) && return true
    end
    return false
end

function set_iterator!(iterator::GridIterator, ri, zi, rf, zf)
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
        #if point_outside_grid(iterator.r_range, iterator.z_range, ri, zi) &&
        #   point_outside_grid(iterator.r_range, iterator.z_range, rf, zf) ||
        #   ((ri == rf) && (zi == zf))
        iterator.ri = ri
        iterator.rf = rf
        iterator.zi = zi
        iterator.zf = zf
        iterator.finished = true
        return
    end
    #println("ri $ri zi $zi rf $rf zf $zf")
    ri, zi = get_intersection_with_grid(
        ri,
        zi,
        iterator.r_range,
        iterator.z_range,
        ri,
        rf,
        zi,
        zf,
        from_outside = true,
    )
    #println("ri $ri zi $zi rf $rf zf $zf")
    rf, zf = get_intersection_with_grid(
        rf,
        zf,
        iterator.r_range,
        iterator.z_range,
        ri,
        rf,
        zi,
        zf,
        from_outside = false,
    )
    direction_r = sign(rf - ri)
    direction_z = sign(zf - zi)
    current_r_idx = searchsorted_first(iterator.r_range, ri, direction_r)
    current_z_idx = searchsorted_first(iterator.z_range, zi, direction_z)
    #println("ri $ri zi $zi rf $rf zf $zf")
    #println("----------------------------")
    #println("current_r $(current_r_idx)")
    #println("current_z $(current_z_idx)")
    iterator.ri = ri
    iterator.rf = rf
    iterator.zi = zi
    iterator.zf = zf
    iterator.direction_r = direction_r
    iterator.direction_z = direction_z
    iterator.current_r_idx = current_r_idx
    iterator.current_z_idx = current_z_idx
    iterator.intersection[1] = ri #iterator.r_range[current_r_idx]
    iterator.intersection[2] = zi #iterator.z_range[current_z_idx]
    iterator.finished = false
    return
end


current_r(iterator::GridIterator) = iterator.r_range[iterator.current_r_idx]
current_z(iterator::GridIterator) = iterator.z_range[iterator.current_z_idx]
next_r(iterator::GridIterator) =
    iterator.r_range[iterator.current_r_idx + iterator.direction_r]
next_z(iterator::GridIterator) =
    iterator.z_range[iterator.current_z_idx + iterator.direction_z]
previous_r(iterator::GridIterator) =
    iterator.r_range[iterator.current_r_idx - iterator.direction_r]
previous_z(iterator::GridIterator) =
    iterator.z_range[iterator.current_z_idx - iterator.direction_z]

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

function step_ray!(iterator::GridIterator, lambda_r, lambda_z)
    if lambda_r < lambda_z || lambda_z === NaN
        iterator.current_r_idx += iterator.direction_r
        lambda = lambda_r
    elseif lambda_z < lambda_r || lambda_r === NaN
        iterator.current_z_idx += iterator.direction_z
        lambda = lambda_z
    else
        iterator.current_r_idx += iterator.direction_r
        iterator.current_z_idx += iterator.direction_z
        lambda = lambda_r
    end
end


function next_intersection!(iterator::GridIterator)
    if iterator.finished
        return nothing
    end
    epsilon = 1e-10
    r0 = iterator.intersection[1] #current_r(iterator)
    z0 = iterator.intersection[2] #current_z(iterator)
    r1 = next_r(iterator)
    z1 = next_z(iterator)
    lambda_r = get_time_to_intersection(r1, iterator.ri, iterator.rf)
    lambda_z = get_time_to_intersection(z1, iterator.zi, iterator.zf)
    #println("---")
    #println("r0 $(r0-epsilon) rf $(iterator.rf) r1 $(r1+epsilon)")
    #println("z0 $(z0-epsilon) zf $(iterator.zf) z1 $(z1+epsilon)")
    if (
        (r0 - epsilon <= iterator.rf <= r1 + epsilon) ||
        (r1 - epsilon <= iterator.rf <= r0 + epsilon) ||
        lambda_r == Inf
    ) && (
        (z0 - epsilon <= iterator.zf <= z1 + epsilon) ||
        (z1 - epsilon <= iterator.zf <= z0 + epsilon) ||
        lambda_z == Inf
    )
        iterator.finished = true
        iterator.intersection[1] = iterator.rf
        iterator.intersection[2] = iterator.zf
        #println("inter $(iterator.intersection)")
        return
    end
    lambda = step_ray!(iterator, lambda_r, lambda_z)
    iterator.intersection[1] = iterator.ri + lambda * (iterator.rf - iterator.ri)
    iterator.intersection[2] = iterator.zi + lambda * (iterator.zf - iterator.zi)
    return
end

function get_density(wi::WindInterpolator, iterator::GridIterator)
    return wi.grid.grid[iterator.current_r_idx, iterator.current_z_idx]
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

function get_density(grid::InterpolationGrid, iterator::GridIterator)
    return grid.grid[iterator.current_r_idx, iterator.current_z_idx]
end

function dist(p1::Vector{Float64}, p2::Vector{Float64}, Rg::Float64)
    return sqrt((p1[1] - p2[1])^2 + (p1[2] - p2[2])^2) * Rg
end

function compute_xray_tau(
    grid::InterpolationGrid{T},
    iterator::GridIterator{T},
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
    set_iterator!(iterator, ri, zi, rf, zf)
    initial_point = [ri, zi]
    previous_point = copy(iterator.intersection)
    ret = 0.0
    while !iterator.finished
        previous_point[1] = iterator.intersection[1]
        previous_point[2] = iterator.intersection[2]
        cell_density = get_density(grid, iterator)
        #cell_density = get_density(grid, iterator.intersection)
        distance_from_source = dist(initial_point, iterator.intersection, Rg)
        next_intersection!(iterator)
        intersection_size = dist(previous_point, iterator.intersection, Rg)
        ret = compute_xray_tau_cell(
            intersection_size,
            distance_from_source,
            cell_density,
            ret,
            xray_luminosity,
            Rg,
        )
        #if ret > 40
        #    return ret
        #end
    end
    return ret
end

function compute_xray_tau(grid::InterpolationGrid, ri, zi, rf, zf, xray_luminosity, Rg)
    compute_xray_tau(
        grid::InterpolationGrid,
        grid.iterator,
        ri,
        zi,
        rf,
        zf,
        xray_luminosity,
        Rg,
    )
end

compute_uv_tau_cell(intersection_size, cell_density)::Float64 =
    intersection_size * cell_density * SIGMA_T

function compute_uv_tau(grid::InterpolationGrid, iterator::CellIterator, ri, zi, rf, zf, Rg)
    if grid.grid === nothing
        return 0.0
    end
    #iterator = grid.iterator
    set_iterator!(iterator, ri, zi, rf, zf)
    previous_point = copy(iterator.intersection)
    ret = 0.0
    while !iterator.finished
        previous_point[1] = iterator.intersection[1]
        previous_point[2] = iterator.intersection[2]
        #cell_density = get_density(grid, iterator.intersection)
        cell_density = get_density(grid, iterator)
        next_intersection!(iterator)
        intersection_size = dist(previous_point, iterator.intersection, Rg)
        ret += compute_uv_tau_cell(intersection_size, cell_density)
    end
    return ret
end

compute_uv_tau(grid::InterpolationGrid, ri, zi, rf, zf, Rg) =
    compute_uv_tau(grid::InterpolationGrid, grid.iterator, ri, zi, rf, zf, Rg)
