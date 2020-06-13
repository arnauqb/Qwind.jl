using Interpolations, Roots
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

## X-Ray phyiscs

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
If we want to solve the equation ξ(r_x) = ξ_0, then we need to solve
the implicit equation ξ_0 = L_x / (n r_x^2) * exp(-σ * n * (r_x - r_in)),
where we introduce the change t = log(r_x) to make it less suitable to roundoff
errors.

# Parameters
- xray_luminosity (in cgs)
- number_density (in cgs)
- r_in initial material radius (in cgs)
- r_x ionization radius candidate (in cgs)
- target ionization parameter (in cgs)
"""
function ionization_radius_kernel(
    xray_luminosity,
    number_density,
    r_in,
    r_x,
    target_ionization_parameter,
)
    return log(
        xray_luminosity /
        (target_ionization_parameter * number_density * r_x^2),
    ) - number_density * SIGMA_T * (r_x - r_in)
end

function ionization_radius(
    xray_luminosity,
    number_density,
    r_in,
    target_ionization_parameter = 1e5;
    atol = 1e-2,
    rtol = 0,
)
    # t = log(r_x) => r_x = e^t
    func(t) = ionization_radius_kernel(
        xray_luminosity,
        number_density,
        r_in,
        exp(t),
        target_ionization_parameter,
    )
    t0 = find_zero(func, (-40, 40), Bisection(), atol = atol, rtol = rtol)
    return exp(t0)
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
    r_1 = 0.0,
    z_1 = 0.0,
)
    d = sqrt((r - r_1)^2 + (z - z_1)^2)
    return number_density * xray_opacity(ionization_parameter) * d
end

function xray_tau(
    r::Vector,
    z::Vector,
    number_density::Vector,
    ionization_parameter::Vector,
)
    ret = 0.0
    r1, z1 = 0.0, 0.0
    for i = 1:length(r)
        ret += xray_tau(
            r[i],
            z[i],
            number_density[i],
            ionization_parameter[i],
            r1,
            z1,
        )
        r1, z1 = r[i], z[i]
    end
    return ret
end

## UV Physics


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
        (fm_interpolation_data["k_interp_x"],),
        fm_interpolation_data["k_interp_y"],
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
        (fm_interpolation_data["eta_interp_x"],),
        fm_interpolation_data["eta_interp_y"],
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
