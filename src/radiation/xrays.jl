"""
Ionization parameter ξ
"""
function compute_ionization_parameter(
    r,
    z,
    vr,
    vz,
    number_density,
    tau_x,
    xray_luminosity,
    Rg;
    zh = 0.0,
)
    d = sqrt(r^2 + z^2)# * Rg
    # relativistic correction to the flux
    beta = max(min(1.0, sqrt(vr^2 + vz^2)), 1e-5)
    gamma = 1.0 / sqrt(1 - beta^2)
    cosθ = (r * vr + (z - zh) * vz) / (d * beta)
    flux_correction = 1.0 / (gamma * (1 + beta * cosθ))^4
    return max(
        xray_luminosity * exp(-tau_x) / (number_density * d^2 * Rg^2) * flux_correction,
        1e-20,
    )
end

compute_ionization_parameter(r, z, number_density, tau_x, xray_luminosity, Rg; zh = 0.0) =
    compute_ionization_parameter(
        r,
        z,
        0.0,
        0.0,
        number_density,
        tau_x,
        xray_luminosity,
        Rg,
        zh = zh,
    )

compute_ionization_parameter(radiation::Radiation, r, z, vr, vz, number_density, tau_x) =
    compute_ionization_parameter(
        r,
        z,
        vr,
        vz,
        number_density,
        tau_x,
        radiation.xray_luminosity,
        radiation.bh.Rg,
        zh = radiation.z_disk,
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
