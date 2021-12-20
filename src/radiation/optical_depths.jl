export compute_optical_depth

function ionization_cell_xi_kernel(;
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
            log(
                xray_luminosity /
                (target_ionization_parameter * cell_density * dx^2),
            ) - cell_density * SIGMA_T * (dx - distance_from_source) - taux0
        return ret
    end
end

function compute_optical_depth_cell(
    density_grid::DensityGrid,
    iterator::GridIterator,
    ::BoostOpacity2;
    r_idx,
    z_idx,
    ri,
    zi,
    phii,
    rf,
    zf,
    phif,
    tau_to_cell,
    source_luminosity,
    Rg,
)
    cell_density = density_grid.grid[r_idx, z_idx]
    distance_to_source = compute_distance_cylindrical(
        ri,
        phii,
        zi,
        iterator.previous_point[1],
        iterator.previous_point[2],
        iterator.previous_point[3],
    ) * Rg
    intersection_size = dist_to_intersection(iterator, iterator.previous_point) * Rg
    f(t) = ionization_cell_xi_kernel(
        xray_luminosity = source_luminosity,
        cell_density = cell_density,
        distance_from_source = distance_to_source,
        taux0 = tau_to_cell,
        dx = t,
    )
    x1 = distance_to_source
    x2 = distance_to_source + intersection_size
    f1 = f(x1)
    f2 = f(x2)
    if f1 > 0 && f2 > 0
        return cell_density * intersection_size * SIGMA_T
    elseif f1 < 0 && f2 < 0
        return cell_density * intersection_size * 100 * SIGMA_T
    else
        dx0 = find_zero(f, (x1, x2), Bisection())
        return cell_density *
               (
                   (dx0 - distance_to_source) +
                   100 * (distance_to_source + intersection_size - dx0)
               ) * SIGMA_T
    end
end


function compute_optical_depth_cell(
    density_grid::DensityGrid,
    iterator::GridIterator,
    ::BoostOpacity;
    r_idx,
    z_idx,
    ri,
    zi,
    phii,
    rf,
    zf,
    phif,
    tau_to_cell,
    source_luminosity,
    Rg,
)
    cell_density = density_grid.grid[r_idx, z_idx]
    dx = sqrt(source_luminosity * exp(-tau_to_cell) / cell_density / 1e5)
    distance_to_source = compute_distance_cylindrical(
        ri,
        phii,
        zi,
        iterator.previous_point[1],
        iterator.previous_point[2],
        iterator.previous_point[3],
    ) * Rg
    intersection_size = dist_to_intersection(iterator, iterator.previous_point) * Rg
    ret = cell_density * SIGMA_T
    if distance_to_source > dx
        ret = 100 * ret * intersection_size
    else
        if distance_to_source + intersection_size > dx
            ret =
                ret * (
                    (dx - distance_to_source) +
                    100 * (distance_to_source + intersection_size - dx)
                )
        else
            ret = ret * intersection_size
        end
    end
    return ret
end


function compute_optical_depth_cell(
    density_grid::DensityGrid,
    iterator::GridIterator,
    ::ThomsonOpacity;
    r_idx,
    z_idx,
    ri,
    zi,
    phii,
    rf,
    zf,
    phif,
    tau_to_cell,
    source_luminosity,
    Rg,
)
    cell_density = density_grid.grid[r_idx, z_idx]
    intersection_size = dist_to_intersection(iterator, iterator.previous_point)
    return intersection_size * cell_density * SIGMA_T * Rg
end

function compute_optical_depth(
    density_grid::DensityGrid,
    iterator,
    opacity::OpacityFlag;
    ri,
    phii,
    zi,
    rf,
    phif,
    zf,
    Rg,
    max_tau = 50.0,
    source_luminosity,
)
    set_iterator!(iterator, ri, phii, zi, rf, phif, zf)
    ret = 0.0
    while !iterator.finished
        iterator.previous_point .= iterator.intersection
        r_idx = max(iterator.current_r_idx, 1)
        z_idx = iterator.current_z_idx
        next_intersection!(iterator)
        ret +=
            compute_optical_depth_cell(
                density_grid,
                iterator,
                opacity,
                tau_to_cell = ret,
                r_idx = r_idx,
                z_idx = z_idx,
                ri = ri,
                zi = zi,
                phii = phii,
                rf = rf,
                zf = zf,
                phif = phif,
                source_luminosity = source_luminosity,
                Rg = Rg,
            ) 
        if ret > max_tau
            return max_tau
        end
    end
    return ret
end
