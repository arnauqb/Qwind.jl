export compute_optical_depth

function compute_optical_depth_cell(
    density_grid::DensityGrid,
    ionization_grid::IonizationGrid,
    ::BoostOpacity;
    mu_electron = 1.17,
    mu_nucleon = 0.61,
    r_idx,
    z_idx,
)
    iterator = density_grid.iterator
    cell_density = density_grid.grid[r_idx, z_idx]
    cell_ion = 1e15
    cell_ion = ionization_grid.grid[
        min(r_idx, ionization_grid.nr),
        min(z_idx, ionization_grid.nz),
    ]
    intersection_size = dist_to_intersection(iterator, iterator.previous_point)
    ret = intersection_size * cell_density * mu_electron
    if cell_ion < 1e5
        return 100 * ret
    else
        return ret
    end
end


function compute_optical_depth_cell(
    density_grid::DensityGrid,
    ionization_grid::IonizationGrid,
    ::ThomsonOpacity;
    mu_electron = 1.17,
    mu_nucleon = 0.61,
    r_idx,
    z_idx,
)
    iterator = density_grid.iterator
    cell_density = density_grid.grid[r_idx, z_idx]
    intersection_size = dist_to_intersection(iterator, iterator.previous_point)
    return intersection_size * cell_density * mu_electron
end

function compute_optical_depth(
    iterator,
    density_grid::DensityGrid,
    ionization_grid::IonizationGrid,
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
            ionization_grid,
            opacity,
            mu_electron = mu_electron,
            mu_nucleon = mu_nucleon,
            r_idx = r_idx,
            z_idx = z_idx,
        )
        if ret * sigmarg > max_tau
            return max_tau
        end
    end
    return ret * sigmarg
end
