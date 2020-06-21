using Interpolations, Cubature
export nt_rel_factors,
    compute_ionization_parameter,
    compute_xray_opacity,
    Interp,
    NoInterp,
    integrate_radiation_force_integrand

abstract type InterpolationType <: Flag end
struct NoInterp <: InterpolationType end
struct Interp <: InterpolationType end

"""
Novikov-Thorne relativistic factors for the AD spectrum.
If the passed radius is smaller than ISCO, it returns 0.

Parameters
----------
radius
    disc radius measured in R_g
spin
    normalized black hole spin (-1, 1)
isco
    innermost stable circular orbit measured in Rg
"""
function nt_rel_factors(radius, spin, isco)
    if radius <= isco
        return 0.0
    end
    yms = sqrt(isco)
    y1 = 2 * cos((acos(spin) - pi) / 3)
    y2 = 2 * cos((acos(spin) + pi) / 3)
    y3 = -2 * cos(acos(spin) / 3)
    y = sqrt(radius)
    C = 1 - 3 / radius + 2 * spin / radius^1.5
    B =
        3 * (y1 - spin)^2 * log((y - y1) / (yms - y1)) /
        (y * y1 * (y1 - y2) * (y1 - y3))
    B +=
        3 * (y2 - spin)^2 * log((y - y2) / (yms - y2)) /
        (y * y2 * (y2 - y1) * (y2 - y3))
    B +=
        3 * (y3 - spin)^2 * log((y - y3) / (yms - y3)) /
        (y * y3 * (y3 - y1) * (y3 - y2))
    A = 1 - yms / y - 3 * spin * log(y / yms) / (2 * y)
    factor = (A - B) / C
    return factor
end

function nt_rel_factors(radiation::Radiation, radius)
    return nt_rel_factors(radius, radiation.bh.spin, radiation.bh.isco)
end
"""
Ionization parameter ξ
"""
function compute_ionization_parameter(
    r,
    z,
    number_density,
    tau_x,
    xray_luminosity,
    Rg,
)
    d = sqrt(r^2 + z^2) * Rg
    return max(xray_luminosity * exp(-tau_x) / (number_density * d^2), 1e-20)
end

compute_ionization_parameter(
    radiation::Radiation,
    r,
    z,
    number_density,
    tau_x,
) = compute_ionization_parameter(
    r,
    z,
    number_density,
    tau_x,
    radiation.xray_luminosity,
    radiation.Rg,
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
compute_force_multiplier_k(ionization_parameter, mode::Interp) =
    force_multiplier_k_log_interpolator(log10(ionization_parameter))
compute_force_multiplier_k(ionization_parameter) =
    compute_force_multiplier_k(ionization_parameter, Interp())
function compute_force_multiplier_k(ionization_parameter, mode::NoInterp)
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
compute_force_multiplier_eta(ionization_parameter, mode::Interp) =
    10 .^ force_multiplier_eta_log_interpolator(log10(ionization_parameter))
compute_force_multiplier_eta(ionization_parameter) =
    compute_force_multiplier_eta(ionization_parameter, Interp())

function compute_force_multiplier_eta(ionization_parameter, mode::NoInterp)
    if (log10(ionization_parameter) < 0.5)
        aux = 6.9 * exp(0.16 * ionization_parameter^0.4)
        eta_max = 10^aux
    else
        aux = 9.1 * exp(-7.96e-3 * ionization_parameter)
        eta_max = 10^aux
    end
end

"This is the sobolev optical depth parameter for the force multiplier"
function compute_tau_eff(number_density, dv_dr)
    v_thermal_sk = 6.77652505049944e-5 # thermal velocity at T=25e3 K in C units
    if dv_dr == 0
        return 1.0
    end
    @assert number_density >= 0
    t = number_density * SIGMA_T * abs(v_thermal_sk / dv_dr)
    return t
end


"""
Computes the analytical approximation for the force multiplier,
from Stevens and Kallman 1990. Note that we modify it slightly to avoid
numerical overflow.
"""
function compute_force_multiplier(
    t,
    ionization_parameter,
    mode::InterpolationType,
)
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
    compute_force_multiplier(t, ionization_parameter, Interp())

"""
Computes the disc radiation field at the point (r,z) by performing
an integral over the disc.

# Parameters
- r: radial coordinate [Rg]
- z: height coordinate [Rg]
- radiative_efficiency: accretion radiative efficiency
"""
function compute_disc_radiation_field(
    radiation::Radiation,
    r,
    z,
    radiative_efficiency,
    r_min,
    r_max = 1600,
    atol = 0,
    rtol = 1e-4,
    norm = Cubature.INDIVIDUAL,
    maxevals = 0,
)
    res, err =
        integrate_radiation_force_integrand(radiation, r, z, r_min, r_max)
    radiation_constant = 3 / (8 * pi * radiation.radiative_efficiency)
    force = radiation_constant .* res
    return force
end

"""
Integrates the radiation acceleration integrand on the disc.

# Parameters
- r : radial coordinate [Rg]
- z : height coordinate [Rg]
-
"""
function integrate_radiation_force_integrand(
    radiation::Radiation,
    r,
    z,
    r_min,
    r_max = 1600;
    phi_min = 0.0,
    phi_max = π,
    atol = 0,
    rtol = 1e-4,
    norm = Cubature.INDIVIDUAL,
    maxevals = 0,
)
    f(x, v) = radiation_force_integrand!(radiation, v, x[1], x[2], r, z)
    return hcubature(
        2,
        f,
        (r_min, phi_min),
        (r_max, phi_max),
        abstol = atol,
        reltol = rtol,
        error_norm = norm,
        maxevals = maxevals,
    )
end
