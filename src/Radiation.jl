using Interpolations, Roots
export compute_xray_opacity,
    compute_xray_tau,
    compute_uv_tau,
    compute_uv_opacity,
    compute_ionization_parameter,
    compute_ionization_radius,
    compute_force_multiplier,
    compute_tau_eff,
    compute_radiation_force,
    xray_luminosity,
    bolometric_luminosity,
    eddington_luminosity,
    mass_accretion_rate,
    Radiation

export NoInterp, Interp, QuadTree, NoGrid
abstract type GridType <: Flag end
abstract type InterpolationType <: Flag end
struct NoInterp <: InterpolationType end
struct Interp <: InterpolationType end
struct QuadTree <: GridType end
struct NoGrid <: GridType end


struct Radiation
    bh::BlackHole
    f_uv::Float64
    f_x::Float64
    uv_fractions::Array{Float64, 1}
end

xray_luminosity(radiation::Radiation) = radiation.f_x * bolometric_luminosity(radiation.bh)
bolometric_luminosity(radiation::Radiation) = bolometric_luminosity(radiation.bh)
eddington_luminosity(radiation::Radiation) = eddington_luminosity(radiation.bh)
mass_accretion_rate(radiation::Radiation) = mass_accretion_rate(radiation.bh)


## X-Ray phyiscs

"""
Ionization parameter ξ
"""
function compute_ionization_parameter(r, z, number_density, tau_x, xray_luminosity)
    d = sqrt(r^2 + z^2)
    return xray_luminosity * exp(-tau_x) / (number_density * d^2)
end

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

function ionization_radius_kernel(
    xray_luminosity,
    number_density,
    r_in,
    r_x,
    target_ionization_parameter,
)
    if r_x < r_in
        return log(xray_luminosity / (target_ionization_parameter * 1e2 * r_x^2))
    else
        return log(xray_luminosity / (target_ionization_parameter * number_density * r_x^2)) -
            number_density * SIGMA_T * (r_x - r_in)
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
- target ionization parameter (in cgs)

# Returns
- r_x ionization radius (in cgs)
"""
function compute_ionization_radius(
    xray_luminosity,
    number_density,
    r_in,
    target_ionization_parameter = 1e5;
    atol = 1e-4,
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
function compute_xray_tau(r, z, number_density, ionization_parameter, r_1 = 0.0, z_1 = 0.0)
    d = sqrt((r - r_1)^2 + (z - z_1)^2)
    return number_density * compute_xray_opacity(ionization_parameter) * d
end

function compute_xray_tau(
    r::Vector,
    z::Vector,
    number_density::Vector,
    ionization_parameter::Vector,
)
    ret = 0.0
    r1, z1 = 0.0, 0.0
    for i = 1:length(r)
        ret +=
            compute_xray_tau(r[i], z[i], number_density[i], ionization_parameter[i], r1, z1)
        r1, z1 = r[i], z[i]
    end
    return ret
end

"""
Risaliti & Elvis (2010) implementation of the X-Ray optical depth.

# Parameters
- r : radial position (cm)
- z : height position (cm)
- shielding density : average atmospheric number density (cm^-3)
- local density : local density at the position (cm^-3)
- r_in : initial radius of the material(wind) (cm)
- r_x : ionization radius (cm)
- r_0 : initial radius of the current streamline (cm)
- mode : RE mode
"""
function compute_xray_tau(r, z, shielding_density, local_density, r_in, r_x, r_0; mode::RE)
    if r <= r_in
        return 0.0
    end
    if r <= r_0
        tau_r = 0.0
        if r_x < r
            tau_r_0 = shielding_density * ((r_x - r_in) + (r - r_x) * 100)
        else
            tau_r_0 = shielding_density * (r - r_in)
        end
    else
        if r_x < r
            if r_0 < r_x
                tau_r_0 = shielding_density * (r_0 - r_in)
                if r < r_x
                    tau_r = local_density * (r - r_x)
                else
                    tau_r = local_density * ((r_x - r_0) + (r - r_x) * 100)
                end
            else
                tau_r_0 = shielding_density * ((r_x - r_in) + 100 * (r_0 - r_x))
                tau_r = local_density * (r - r_0) * 100
            end
        else
            tau_r_0 = shielding_density * (r_0 - r_in)
            tau_r = local_density * (r - r_0)
        end
    end
    upper_projection = 1.0 / cos(r / sqrt(r^2 + z^2))
    return (tau_r + tau_r_0) * SIGMA_T * upper_projection
end


function compute_xray_tau(r, z, density, p::Parameters, mode::RE)
    return compute_xray_tau(
        r,
        z,
        number_density_0(p.line),
        density,
        p.r_in,
        p.r_x,
        r_0(p.line),
        RE(),
    )
end


## UV Physics


"""
UV opacity
"""
compute_uv_opacity() = SIGMA_T

function compute_uv_tau(r, z, number_density, r_0 = 0.0, z_0 = 0.0)
    d = sqrt((r - r_0)^2 + (z - z_0)^2)
    return number_density * compute_uv_opacity() * d
end

"""
Risaliti & Elvis (2010) implementation of the X-Ray optical depth.

# Parameters
- r : radial position (cm)
- z : height position (cm)
- shielding density : average atmospheric number density (cm^-3)
- local density : local density at the position (cm^-3)
- r_in : initial radius of the material(wind) (cm)
- r_x : ionization radius (cm)
- r_0 : initial radius of the current streamline (cm)
- mode : RE mode
"""
function compute_uv_tau(r, z, shielding_density, local_density, r_in, r_0; mode::RE)
    if r <= r_in
        return 0.0
    end
    if r <= r_0
        tau_r = 0.0
        tau_r_0 = shielding_density * (r - r_in)
    else
        tau_r_0 = shielding_density * (r_0 - r_in)
        tau_r = local_density * (r - r_0)
    end
    upper_projection = 1.0 / cos(r / sqrt(r^2 + z^2))
    return (tau_r + tau_r_0) * SIGMA_T * upper_projection
end

function compute_uv_tau(r, z, density, p::Parameters, mode::RE)
    return compute_uv_tau(
        r,
        z,
        number_density_0(p.line),
        density,
        p.r_in,
        r_0(p.line),
        RE(),
    )
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
function compute_tau_eff(number_density, dv_dr, v_th)
    if dv_dr == 0
        return 1.
    end
    @assert density >= 0
    t = number_density * SIGMA_T * abs(v_th / dv_dr)
    return t
end


"""
Computes the analytical approximation for the force multiplier,
from Stevens and Kallman 1990. Note that we modify it slightly to avoid
numerical overflow.
"""
function compute_force_multiplier(t, ionization_parameter, mode::InterpolationType)
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

function compute_radiation_force(
    integrand,
    r,
    z,
    force_multiplier,
    M,
    r_lims = (6.0, 1500.0),
    phi_lims = (0.0, 2π),
)
    res, err = integrate_2d(integrand, x_lims = r_lims, y_lims = phi_lims)
    radiation_constant = 3 * G * M / (8π)
    return (1 + force_multiplier) * radiation_constant * res
end
