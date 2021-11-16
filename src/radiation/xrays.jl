export compute_ionization_parameter, compute_xray_opacity

"""
Ionization parameter ξ
"""
function compute_ionization_parameter(;
    r,
    z,
    vr,
    vz,
    number_density,
    tau_x,
    xray_luminosity,
    Rg,
    include_scattering = false,
    scattered_luminosity_grid = nothing,
    density_grid = nothing,
    absorption_opacity = Boost(),
    zh = 0.0,
    mu_electron = 1.17,
    mu_nucleon = 0.61,
)
    d = sqrt(r^2 + z^2)# * Rg
    # relativistic correction to the flux
    beta = max(min(1.0, sqrt(vr^2 + vz^2)), 1e-5)
    gamma = 1.0 / sqrt(1 - beta^2)
    cosθ = (r * vr + (z - zh) * vz) / (d * beta)
    flux_correction = 1.0 / (gamma * (1 + beta * cosθ))^4
    ret = xray_luminosity * exp(-tau_x) / (number_density * d^2 * Rg^2) * flux_correction
    if include_scattering
        ret +=
            scattered_flux_in_point(
                scattered_luminosity_grid,
                density_grid,
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

compute_ionization_parameter(radiation::Radiation; r, z, vr, vz, number_density, tau_x) =
    compute_ionization_parameter(
        r = r,
        z = z,
        vr = vr,
        vz = vz,
        number_density = number_density,
        tau_x = tau_x,
        xray_luminosity = radiation.xray_luminosity,
        Rg = radiation.bh.Rg,
        include_scattering = radiation.xray_scattering,
        density_grid = radiation.wi.density_grid,
        absorption_opacity = radiation.xray_opacity,
        zh = radiation.z_xray,
        mu_electron = radiation.mu_electron,
        mu_nucleon = radiation.mu_nucleon,
        scattered_luminosity_grid = radiation.wi.scattered_lumin_grid,
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

