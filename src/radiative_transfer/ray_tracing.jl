using Distances
export compute_cell_intersection, next!, GridIterator, compute_xray_tau, compute_uv_tau

function point_outside_grid(r_range, z_range, r, z)
    if !(r_range[1] <= r <= r_range[end]) || !(z_range[1] <= z <= z_range[end])
        return true
    else
        return false
    end
end

get_time_to_intersection(x, xi, xf) = abs((x - xi) / (xf - xi))

function get_intersection_with_grid(r, z, r_range, z_range, ri, rf, zi, zf; from_outside)
    lambda_r = Inf
    lambda_z = Inf
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
        lambda = max(lambda_r, lambda_z)
    else
        lambda = min(lambda_r, lambda_z)
    end
    if lambda != Inf && !isnan(lambda)
        r = ri + lambda * (rf - ri)
        z = zi + lambda * (zf - zi)
    end
    return r, z
end

mutable struct GridIterator <: CellIterator
    r_range::Vector{Float64}
    z_range::Vector{Float64}
    ri::Float64
    zi::Float64
    rf::Float64
    zf::Float64
    direction_r::Int
    direction_z::Int
    current_r_idx::Int
    current_z_idx::Int
    state::Int
    finished::Bool
    function GridIterator(r_range, z_range, ri, zi, rf, zf)
        if point_outside_grid(r_range, z_range, ri, zi) &&
           point_outside_grid(r_range, z_range, rf, zf)
            return new(r_range, z_range, ri, zf, rf, zf, 1, 1, 1, 1, 1, true)
        end
        #println("ri $ri zi $zi rf $rf zf $zf")
        ri, zi = get_intersection_with_grid(
            ri,
            zi,
            r_range,
            z_range,
            ri,
            rf,
            zi,
            zf,
            from_outside = true,
        )
        rf, zf = get_intersection_with_grid(
            rf,
            zf,
            r_range,
            z_range,
            ri,
            rf,
            zi,
            zf,
            from_outside = false,
        )
        direction_r = sign(rf - ri)
        direction_z = sign(zf - zi)
        current_r_idx = searchsorted_first(r_range, ri, direction_r)
        current_z_idx = searchsorted_first(z_range, zi, direction_z)
        #println("ri $ri zi $zi rf $rf zf $zf")
        #println("----------------------------")
        return new(
            r_range,
            z_range,
            ri,
            zi,
            rf,
            zf,
            direction_r,
            direction_z,
            current_r_idx,
            current_z_idx,
            0,
            false,
        )
    end
end

current_r(iterator::GridIterator) = iterator.r_range[iterator.current_r_idx]
current_z(iterator::GridIterator) = iterator.z_range[iterator.current_z_idx]
next_r(iterator::GridIterator) =
    iterator.r_range[iterator.current_r_idx + iterator.direction_r]
next_z(iterator::GridIterator) =
    iterator.z_range[iterator.current_z_idx + iterator.direction_z]

function GridIterator(interpolator::DensityInterpolator, ri, zi, rf, zf)
    return GridIterator(interpolator.grid.r_range, interpolator.grid.z_range, ri, zi, rf, zf)
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


function next_intersection!(iterator::GridIterator; init = false)
    if iterator.finished
        return nothing
    end
    if init
        return [iterator.ri, iterator.zi], iterator.state
    end
    epsilon = 1e-10
    r0 = current_r(iterator)
    z0 = current_z(iterator)
    r1 = next_r(iterator)
    z1 = next_z(iterator)
    lambda_r = get_time_to_intersection(r1, iterator.ri, iterator.rf)
    lambda_z = get_time_to_intersection(z1, iterator.zi, iterator.zf)
    #println("r0 $(r0-epsilon) rf $(iterator.rf) r1 $(r1+epsilon)")
    #println("z0 $(z0-epsilon) rf $(iterator.zf) z1 $(z1+epsilon)")
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
        return [iterator.rf, iterator.zf], iterator.state
    end
    lambda = step_ray!(iterator, lambda_r, lambda_z)
    iterator.state += 1
    intersection = [
        iterator.ri + lambda * (iterator.rf - iterator.ri),
        iterator.zi + lambda * (iterator.zf - iterator.zi),
    ]
    return (intersection, iterator.state)
end

Base.iterate(iterator::GridIterator) = next_intersection!(iterator, init = true)
Base.iterate(iterator::GridIterator, state) = next_intersection!(iterator)
Base.eltype(::Type{GridIterator}) = Array{Float64,1}

function compute_xray_tau_cell(
    intersection_size,
    distance_from_source,
    cell_density,
    taux0,
    xray_luminosity,
)
    column_density = cell_density * intersection_size
    xi0 = xray_luminosity / (cell_density * distance_from_source^2)
    f(t) = t - log(xi0) + taux0 + min(40, column_density * compute_xray_opacity(exp(t)))
    if f(20) < 0
        xi = xi0
    elseif f(-20) > 0
        xi = 1e-20
    else
        t = find_zero(f, (-20, 20), Bisection(), atol = 0, rtol = 0.1)
        xi = exp(t)
    end
    taux = compute_xray_opacity(xi) * column_density
    return taux
end

function compute_xray_tau(
    density_interpolator::DensityInterpolator,
    ri,
    zi,
    rf,
    zf,
    xray_luminosity,
    Rg,
)
    iterator = GridIterator(
        density_interpolator.grid.r_range,
        density_interpolator.grid.z_range,
        ri,
        zi,
        rf,
        zf,
    )
    initial_point = [ri, zi]
    previous_point = [current_r(iterator), current_z(iterator)]
    ret = 0
    for current_point in iterator
        distance_from_source = Distances.evaluate(Euclidean(), initial_point, current_point) * Rg
        intersection_size = Distances.evaluate(Euclidean(), previous_point, current_point) * Rg
        cell_density =
            get_density(density_interpolator, previous_point[1], previous_point[2])
        previous_point = current_point
        ret += compute_xray_tau_cell(
            intersection_size,
            distance_from_source,
            cell_density,
            ret,
            xray_luminosity,
        )
        if ret > 40
            return ret
        end
    end
    return ret
end

compute_uv_tau_cell(intersection_size, cell_density) =
    intersection_size * cell_density * SIGMA_T

function compute_uv_tau(density_interpolator::DensityInterpolator, ri, zi, rf, zf, Rg)
    iterator = GridIterator(
        density_interpolator.grid.r_range,
        density_interpolator.grid.z_range,
        ri,
        zi,
        rf,
        zf,
    )
    previous_point = [current_r(iterator), current_z(iterator)]
    ret = 0
    for current_point in iterator
        intersection_size = Distances.evaluate(Euclidean(), previous_point, current_point) * Rg
        cell_density =
            get_density(density_interpolator, previous_point[1], previous_point[2])
        previous_point = current_point
        ret += compute_uv_tau_cell(intersection_size, cell_density)
        if ret > 20
            return ret
        end
    end
    return ret
end