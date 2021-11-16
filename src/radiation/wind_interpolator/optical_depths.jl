export compute_optical_depth

function compute_optical_depth_cell(
    density_grid::DensityGrid,
    #ionization_grid::IonizationGrid,
    ::BoostOpacity;
    mu_electron = 1.17,
    mu_nucleon = 0.61,
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
    iterator = density_grid.iterator
    cell_density = density_grid.grid[r_idx, z_idx]
    dx = sqrt(source_luminosity * exp(-tau_to_cell) / cell_density / 1e5) / Rg
    distance_to_source = compute_distance_cylindrical(
        ri,
        phii,
        zi,
        iterator.previous_point[1],
        iterator.previous_point[2],
        iterator.previous_point[3],
    )
    intersection_size = dist_to_intersection(iterator, iterator.previous_point)
    ret = cell_density * mu_electron
    if distance_to_source > dx
        ret = 100 * ret * intersection_size
    else
        if distance_to_source + intersection_size > dx
            ret = ret * (dx - distance_to_source) +
                   100 * (distance_to_source + intersection_size - dx)
        else
            ret = ret * intersection_size
        end
    end
    return ret
    #delta_boost = max(distance_to_source - dx, 0)
    #ret = (delta_boost * 100 + dx) * mu_electron * cell_density
    #return ret
    #ret = intersection_size * cell_density * mu_electron
    #if cell_ion < 1e5
    #    return 100 * ret
    #else
    #    return ret
    #end
end


function compute_optical_depth_cell(
    density_grid::DensityGrid,
    #ionization_grid::IonizationGrid,
    ::ThomsonOpacity;
    mu_electron = 1.17,
    mu_nucleon = 0.61,
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
    iterator = density_grid.iterator
    cell_density = density_grid.grid[r_idx, z_idx]
    intersection_size = dist_to_intersection(iterator, iterator.previous_point)
    return intersection_size * cell_density * mu_electron
end

function compute_optical_depth(
    iterator,
    density_grid::DensityGrid,
    #ionization_grid::IonizationGrid,
    opacity::OpacityFlag;
    ri,
    phii,
    zi,
    rf,
    phif,
    zf,
    Rg,
    mu_electron = 1.17,
    mu_nucleon = 0.61,
    max_tau = 50,
    source_luminosity,
)
    set_iterator!(iterator, ri, phii, zi, rf, phif, zf)
    ret = 0.0
    sigmarg = Rg * SIGMA_T
    while !iterator.finished
        iterator.previous_point .= iterator.intersection
        r_idx = max(iterator.current_r_idx, 1)
        z_idx = iterator.current_z_idx
        next_intersection!(iterator)
        ret += compute_optical_depth_cell(
            density_grid,
            #ionization_grid,
            opacity,
            mu_electron = mu_electron,
            mu_nucleon = mu_nucleon,
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
        ) * sigmarg
        if ret > max_tau
            return max_tau
        end
    end
    return ret 
end
