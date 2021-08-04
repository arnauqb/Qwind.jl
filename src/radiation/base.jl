using Interpolations, Cubature
export disk_nt_rel_factors,
    compute_force_multiplier,
    compute_ionization_parameter,
    compute_xray_opacity,
    integrate_radiation_force_integrand



function disk_nt_rel_factors(radiation::Radiation, radius)
    return disk_nt_rel_factors(radius, radiation.spin, radiation.isco)
end

function compute_eddington_luminosity(radiation::Radiation)
    return 4 * π * radiation.Rg * C^3 * M_P / SIGMA_T
end

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
    d = sqrt(r^2 + z^2) * Rg
    # relativistic correction to the flux
    beta = max(sqrt(vr^2 + vz^2), 1e-6)
    gamma = 1.0 / sqrt(1 - beta^2)
    cosθ = (r * vr + (z - zh) * vz) / (d * beta)
    flux_correction = 1.0 / (gamma * (1 + beta * cosθ))^4
    return max(
        xray_luminosity * exp(-tau_x) / (number_density * d^2) * flux_correction,
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

compute_ionization_parameter(radiation::Radiation, r, z, number_density, tau_x) =
    compute_ionization_parameter(
        r,
        z,
        number_density,
        tau_x,
        radiation.xray_luminosity,
        radiation.Rg,
        zh = radiation.zh,
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


"UV opacity for optical depths calculations is assumed to just be σ"
compute_uv_opacity() = SIGMA_T

# Force multiplier

"""
Calculates the fitting parameter k for the force multiplier calculation
in Steven & Kallman 1990
"""
force_multiplier_k_log_interpolator = extrapolate(
    interpolate(
        (fm_interpolation_data["k_interp_x"],),
        fm_interpolation_data["k_interp_y"],
        Gridded(Linear()),
    ),
    Flat(),
)
compute_force_multiplier_k(ionization_parameter, mode::FMInterp) =
    force_multiplier_k_log_interpolator(log10(ionization_parameter))
compute_force_multiplier_k(ionization_parameter) =
    compute_force_multiplier_k(ionization_parameter, FMInterp())
function compute_force_multiplier_k(ionization_parameter, mode::FMNoInterp)
    k = 0.03 + 0.385 * exp(-1.4 * ionization_parameter^0.6)
    return k
end

"""
Calculates the fitting parameter eta_max for the force multiplier calculation
in Steven & Kallman 1990
"""
force_multiplier_eta_log_interpolator = extrapolate(
    interpolate(
        (fm_interpolation_data["eta_interp_x"],),
        fm_interpolation_data["eta_interp_y"],
        Gridded(Linear()),
    ),
    Flat(),
)
compute_force_multiplier_eta(ionization_parameter, mode::FMInterp) =
    10 .^ force_multiplier_eta_log_interpolator(log10(ionization_parameter))
compute_force_multiplier_eta(ionization_parameter) =
    compute_force_multiplier_eta(ionization_parameter, FMInterp())

function compute_force_multiplier_eta(ionization_parameter, mode::FMNoInterp)
    if (log10(ionization_parameter) < 0.5)
        aux = 6.9 * exp(0.16 * ionization_parameter^0.4)
        eta_max = 10^aux
    else
        aux = 9.1 * exp(-7.96e-3 * ionization_parameter)
        eta_max = 10^aux
    end
end

"This is the sobolev optical depth parameter for the force multiplier"
function compute_tau_eff(number_density, dv_dr, mu = 0.61, mu_e = 1.17)
    v_thermal_sk = 6.77652505049944e-5 # thermal velocity at T=25e3 K in C units
    if dv_dr == 0
        return 1.0
    end
    @assert number_density >= 0
    t = number_density * SIGMA_T * mu / mu_e * abs(v_thermal_sk / dv_dr)
    return t
end


"""
Computes the analytical approximation for the force multiplier,
from Stevens and Kallman 1990. Note that we modify it slightly to avoid
numerical overflow.
"""
function compute_force_multiplier(t, ionization_parameter, mode::FMInterpolationType)
    @assert t >= 0
    @assert ionization_parameter >= 0
    ALPHA = 0.6
    TAU_MAX_TOL = 1e-3
    k = compute_force_multiplier_k(ionization_parameter, mode)
    eta = compute_force_multiplier_eta(ionization_parameter, mode)
    tau_max = t * eta
    if tau_max < TAU_MAX_TOL
        aux = (1 - ALPHA) * (tau_max^ALPHA)
    else
        aux = ((1 + tau_max)^(1 - ALPHA) - 1) / ((tau_max)^(1 - ALPHA))
    end
    fm = k * t^(-ALPHA) * aux
    @assert fm >= 0
    return fm
end

compute_force_multiplier(t, ionization_parameter) =
    compute_force_multiplier(t, ionization_parameter, FMInterp())

