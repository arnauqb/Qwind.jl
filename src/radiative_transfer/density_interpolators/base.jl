export get_density, update_density_interpolator

get_density(grid::DensityInterpolator, r, z) =
    error("Density Interpolator should implement this function.")

update_density_interpolator(grid::DensityInterpolator, integrators) =
    error("Density Interpolator should implement this function.")

