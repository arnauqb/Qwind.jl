using Interpolations
export xray_opacity,
    xray_tau, uv_tau, uv_opacity, ionization_parameter, force_multiplier
export NoInterp, Interp, QuadTree, NoGrid
abstract type Flag end
abstract type GridType <: Flag end
abstract type InterpolationType <: Flag end
struct NoInterp <: InterpolationType end
struct Interp <: InterpolationType end
struct QuadTree <: GridType end
struct NoGrid <: GridType end

"""
Ionization parameter ξ
"""
function ionization_parameter(r, z, number_density, tau_x, xray_luminosity)
    d = sqrt(r^2 + z^2)
    return xray_luminosity * exp(-tau_x) / (number_density * d^2)
end

"""
X-Ray opacity as a function of ionization parameter.
"""
function xray_opacity(ionization_parameter)
    if ionization_parameter < 1e5
        return 100 * SIGMA_T
    else
        return 1 * SIGMA_T
    end
end

"""
X-Ray optical depth

# Parameters
- r : disc radius (in Rg)
- z : disc height (in Rg)
- number_density: (in cgs)
- ionization_parameter: (in cgs)

"""
function xray_tau(
    r,
    z,
    number_density,
    ionization_parameter,
    r_0 = 0.0,
    z_0 = 0.0,
)
    d = sqrt((r - r_0)^2 + (z - z_0)^2)
    return number_density * xray_opacity(ionization_parameter) * d
end

"""
UV opacity
"""
uv_opacity() = SIGMA_T

function uv_tau(r, z, number_density, r_0 = 0.0, z_0 = 0.0)
    d = sqrt((r - r_0)^2 + (z - z_0)^2)
    return number_density * uv_opacity() * d
end


"""
Calculates the fitting parameter k for the force multiplier calculation
in Steven & Kallman 1990
"""
force_multiplier_k_log_interpolator = extrapolate(
    interpolate(
        (Tables.k_interp_log_xi,),
        Tables.k_interp_k,
        Gridded(Linear()),
    ),
    Flat(),
)
force_multiplier_k(ionization_parameter, mode::Interp) =
    force_multiplier_k_log_interpolator(log10(ionization_parameter))
force_multiplier_k(ionization_parameter) =
    force_multiplier_k(ionization_parameter, Interp())
function force_multiplier_k(ionization_parameter, mode::NoInterp)
    k = 0.03 + 0.385 * exp(-1.4 * ionization_parameter^0.6)
    return k
end

"""
Calculates the fitting parameter eta_max for the force multiplier calculation
in Steven & Kallman 1990
"""
force_multiplier_eta_log_interpolator = extrapolate(
    interpolate(
        (Tables.etamax_interp_log_xi,),
        Tables.etamax_interp_etamax,
        Gridded(Linear()),
    ),
    Flat(),
)
force_multiplier_eta(ionization_parameter, mode::Interp) =
    10 .^ force_multiplier_eta_log_interpolator(log10(ionization_parameter))
force_multiplier_eta(ionization_parameter) =
    force_multiplier_eta(ionization_parameter, Interp())

function force_multiplier_eta(ionization_parameter, mode::NoInterp)
    if (log10(ionization_parameter) < 0.5)
        aux = 6.9 * exp(0.16 * ionization_parameter^0.4)
        eta_max = 10^aux
    else
        aux = 9.1 * exp(-7.96e-3 * ionization_parameter)
        eta_max = 10^aux
    end
end

"""
Computes the analytical approximation for the force multiplier,
from Stevens and Kallman 1990. Note that we modify it slightly to avoid
numerical overflow.
"""
function force_multiplier(t, ionization_parameter, mode::InterpolationType)
    @assert t >= 0
    @assert ionization_parameter >= 0
    ALPHA = 0.6
    TAU_MAX_TOL = 1e-3
    k = force_multiplier_k(ionization_parameter, mode)
    eta = force_multiplier_eta(ionization_parameter, mode)
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

force_multiplier(t, ionization_parameter) =
    force_multiplier(t, ionization_parameter, Interp())
