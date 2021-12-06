export compute_ionization_parameter, compute_xray_opacity

"""
Ionization parameter ξ
"""
function compute_ionization_parameter(
    scattered_luminosity_grid::ScatteredLuminosityGrid,
    iterator::GridIterator,
    density_grid::DensityGrid;
    r,
    z,
    number_density,
    tau_x,
    xray_luminosity,
    Rg,
    scattering_flag,
    absorption_opacity = BoostOpacity(),
    zh = 0.0,
    mu_electron = 1.17,
    mu_nucleon = 0.61,
)
    d = sqrt(r^2 + z^2)# * Rg
    ret = xray_luminosity * exp(-tau_x) / (number_density * d^2 * Rg^2)
    if scattering_flag == Scattering()
        ret +=
            scattered_flux_in_point(
                scattered_luminosity_grid,
                density_grid,
                iterator,
                r = r,
                z = z,
                mu_electron = mu_electron,
                mu_nucleon = mu_nucleon,
                absorption_opacity = absorption_opacity,
                Rg = Rg,
            ) / number_density
    end
    return max(ret, 1e-20)
end

compute_ionization_parameter(
    radiation::Radiation,
    wind::Wind,
    parameters::Parameters;
    r,
    z,
    number_density,
    tau_x,
) = compute_ionization_parameter(
    radiation.scattered_lumin_grid,
    wind.grid_iterator,
    wind.density_grid,
    r = r,
    z = z,
    number_density = number_density,
    tau_x = tau_x,
    xray_luminosity = radiation.xray_luminosity,
    Rg = radiation.bh.Rg,
    scattering_flag = parameters.xray_scattering_flag,
    absorption_opacity = parameters.xray_opacity_flag,
    zh = parameters.z_xray,
    mu_electron = parameters.mu_electron,
    mu_nucleon = parameters.mu_nucleon,
)

"""
X-Ray opacity as a function of ionization parameter.
"""
function compute_xray_opacity(ionization_parameter)
    if ionization_parameter < 1e5
        return 100 * SIGMA_T
    else
        return 1 * SIGMA_T
    end
end

