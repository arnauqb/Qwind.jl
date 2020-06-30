"""
This radiation module implements the radiation quantities as computed
in the first version of Qwind (Risaliti & Elvis 2010, Quera-Bofarull et al. 2020)
"""

using Roots
export RERadiation,
    compute_xray_tau,
    compute_uv_tau,
    radiation_force_integrand!,
    compute_disc_radiation_field

struct RERadiation <: Radiation
    f_uv::Float64
    xray_luminosity::Float64
    shielding_density::Float64
    r_x::Float64
    r_in::Float64
    bh::BlackHole
    function RERadiation(bh::BlackHole, f_uv, f_x, shielding_density, r_in, mu = 0.5)
        bolometric_luminosity = compute_bolometric_luminosity(bh)
        xray_luminosity = f_x * bolometric_luminosity
        r_x = compute_ionization_radius(xray_luminosity, shielding_density, r_in, bh.Rg)
        new(f_uv, xray_luminosity, shielding_density, r_x, r_in, bh)
    end
end

function ionization_radius_kernel(
    xray_luminosity,
    number_density,
    r_in,
    r_x,
    target_ionization_parameter,
    Rg,
)
    r_x = r_x * Rg
    r_in = r_in * Rg
    if r_x < r_in
        return log(xray_luminosity / (target_ionization_parameter * 1e2 * r_x^2))
    else
        return log(
            xray_luminosity / (target_ionization_parameter * number_density * r_x^2),
        ) - number_density * SIGMA_T * (r_x - r_in)
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
    Rg,
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
        Rg,
    )
    t0 = find_zero(func, (-40, 40), Bisection(), atol = atol, rtol = rtol)
    return exp(t0)
end
function compute_ionization_radius(
    radiation::RERadiation,
    number_density,
    r_in,
    target_ionization_parameter = 1e5;
    atol = 1e-4,
    rtol = 0,
)
    return compute_ionization_radius(
        radiation.xray_luminosity,
        number_density,
        r_in,
        radiation.bh.Rg,
        target_ionization_parameter,
        atol,
        rtol,
    )
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
"""
function compute_xray_tau(r, z, shielding_density, local_density, r_in, r_x, r_0, Rg)
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
    return (tau_r + tau_r_0) * SIGMA_T * upper_projection * Rg
end


function compute_xray_tau(radiation::RERadiation, r, z, local_density, r_0)
    return compute_xray_tau(
        r,
        z,
        radiation.shielding_density,
        local_density,
        radiation.r_in,
        radiation.r_x,
        r_0,
        radiation.bh.Rg,
    )
end


"""
Risaliti & Elvis (2010) implementation of the UV optical depth.

# Parameters
- r : radial position (cm)
- z : height position (cm)
- shielding density : average atmospheric number density (cm^-3)
- local density : local density at the position (cm^-3)
- r_in : initial radius of the material(wind) (cm)
- r_0 : initial radius of the current streamline (cm)
"""
function compute_uv_tau(r, z, shielding_density, local_density, r_in, r_0, Rg)
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
    return (tau_r + tau_r_0) * SIGMA_T * upper_projection * Rg
end

function compute_uv_tau(radiation::RERadiation, r, z, local_density, r_0)
    return compute_uv_tau(
        r,
        z,
        radiation.shielding_density,
        local_density,
        radiation.r_in,
        r_0,
        radiation.bh.Rg,
    )
end

"""
Radiation force integrand for the RE model.
"""
function radiation_force_integrand!(radiation::RERadiation, v, r_d, phi_d, r, z)
    nt = nt_rel_factors(radiation, r_d)
    r_projection = (r - r_d * cos(phi_d))
    delta_sq = (r^2 + r_d^2 + z^2 - 2 * r * r_d * cos(phi_d))
    common_projection = 1.0 / (r_d^2 * delta_sq^2)
    v[:] = nt * common_projection * [r_projection, z]
end

"""
Computes the disc radiation field at the point (r,z) by performing
an integral over the disc.

# Parameters
- r: radial coordinate [Rg]
- z: height coordinate [Rg]
- radiative_efficiency: accretion radiative efficiency
"""
function compute_disc_radiation_field(
    radiation::RERadiation,
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
    #println("r : $r, z : $z")
    res, err = integrate_radiation_force_integrand(radiation, r, z, r_min, r_max)
    radiation_constant =
        2 * 3 * radiation.bh.mdot * radiation.f_uv * z / (8 * π * radiation.bh.efficiency)
    force = radiation_constant .* res
    #println("force $force")
    return force
end
