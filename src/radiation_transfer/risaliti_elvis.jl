export RERadiativeTransfer, compute_ionization_radius, compute_xray_tau, compute_uv_tau, compute_disc_radiation_field, radiation_force_integrand!

struct RERadiativeTransfer <: RadiativeTransfer
    radiation::Radiation
    rx::Float64
    rin::Float64
    shielding_density::Float64
    Rg::Float64
    function RERadiativeTransfer(radiation, shielding_density, rin, Rg)
        rx = compute_ionization_radius(radiation, shielding_density, rin, Rg)
        new(radiation, rx, rin, shielding_density, Rg)
    end
end

function ionization_radius_kernel(
    xray_luminosity,
    shielding_density,
    rin,
    rx,
    target_ionization_parameter,
    Rg,
)
    rx = rx * Rg
    rin = rin * Rg
    if rx < rin
        return log(xray_luminosity / (target_ionization_parameter * 1e2 * rx^2))
    else
        return log(
            xray_luminosity / (target_ionization_parameter * shielding_density * rx^2),
        ) - shielding_density * SIGMA_T * (rx - rin)
    end
end

"""
If we want to solve the equation ξ(rx) = ξ_0, then we need to solve
the implicit equation ξ_0 = L_x / (n rx^2) * exp(-σ * n * (rx - rin)),
where we introduce the change t = log(rx) to make it less suitable to roundoff
errors.

# Parameters
- xray_luminosity (in cgs)
- number_density (in cgs)
- rin initial material radius (in cgs)
- target ionization parameter (in cgs)

# Returns
- rx ionization radius (in cgs)
"""
function compute_ionization_radius(
    xray_luminosity,
    shielding_density,
    rin,
    Rg,
    target_ionization_parameter = 1e5;
    atol = 1e-4,
    rtol = 0,
)
    # t = log(rx) => rx = e^t
    func(t) = ionization_radius_kernel(
        xray_luminosity,
        shielding_density,
        rin,
        exp(t),
        target_ionization_parameter,
        Rg,
    )
    t0 = find_zero(func, (-40, 40), Bisection(), atol = atol, rtol = rtol)
    return exp(t0)
end
function compute_ionization_radius(
    radiation::RERadiation,
    shielding_density,
    rin,
    Rg,
    target_ionization_parameter = 1e5;
    atol = 1e-4,
    rtol = 0.0,
)
    return compute_ionization_radius(
        radiation.xray_luminosity,
        shielding_density,
        rin,
        Rg,
        target_ionization_parameter,
        atol = atol,
        rtol = rtol,
    )
end

"""
Risaliti & Elvis (2010) implementation of the X-Ray optical depth.

# Parameters
- r : radial position (cm)
- z : height position (cm)
- shielding density : average atmospheric number density (cm^-3)
- local density : local density at the position (cm^-3)
- rin : initial radius of the material(wind) (cm)
- rx : ionization radius (cm)
- r0 : initial radius of the current streamline (cm)
"""
function compute_xray_tau(r, z, shielding_density, local_density, rin, rx, r0, Rg)
    if r <= rin
        return 0.0
    end
    if r <= r0
        tau_r = 0.0
        if rx < r
            tau_r0 = shielding_density * ((rx - rin) + (r - rx) * 100)
        else
            tau_r0 = shielding_density * (r - rin)
        end
    else
        if rx < r
            if r0 < rx
                tau_r0 = shielding_density * (r0 - rin)
                if r < rx
                    tau_r = local_density * (r - rx)
                else
                    tau_r = local_density * ((rx - r0) + (r - rx) * 100)
                end
            else
                tau_r0 = shielding_density * ((rx - rin) + 100 * (r0 - rx))
                tau_r = local_density * (r - r0) * 100
            end
        else
            tau_r0 = shielding_density * (r0 - rin)
            tau_r = local_density * (r - r0)
        end
    end
    upper_projection = 1.0 / cos(r / sqrt(r^2 + z^2))
    return (tau_r + tau_r0) * SIGMA_T * upper_projection * Rg
end


function compute_xray_tau(radiative_transfer::RERadiativeTransfer, r, z, local_density, r0)
    return compute_xray_tau(
        r,
        z,
        radiative_transfer.shielding_density,
        local_density,
        radiative_transfer.rin,
        radiative_transfer.rx,
        r0,
        radiative_transfer.Rg,
    )
end


"""
Risaliti & Elvis (2010) implementation of the UV optical depth.

# Parameters
- r : radial position (cm)
- z : height position (cm)
- shielding density : average atmospheric number density (cm^-3)
- local density : local density at the position (cm^-3)
- rin : initial radius of the material(wind) (cm)
- r0 : initial radius of the current streamline (cm)
"""
function compute_uv_tau(r, z, shielding_density, local_density, rin, r0, Rg)
    if r <= rin
        return 0.0
    end
    if r <= r0
        tau_r0 = shielding_density * (r - rin)
        tau_r = 0.0
    else
        tau_r0 = shielding_density * (r0 - rin)
        tau_r = local_density * (r - r0)
    end
    upper_projection = 1.0 / cos(r / sqrt(r^2 + z^2))
    return (tau_r + tau_r0) * SIGMA_T * upper_projection * Rg
end

function compute_uv_tau(radiative_transfer::RERadiativeTransfer, r, z, local_density, r0)
    return compute_uv_tau(
        r,
        z,
        radiative_transfer.shielding_density,
        local_density,
        radiative_transfer.rin,
        r0,
        radiative_transfer.Rg,
    )
end

"""
Radiation force integrand for the RE model.
"""
function radiation_force_integrand!(
    radiative_transfer::RERadiativeTransfer,
    v,
    rd,
    phid,
    r,
    z,
)
    nt = nt_rel_factors(radiative_transfer.radiation, rd)
    r_projection = (r - rd * cos(phid))
    delta_sq = (r^2 + rd^2 + z^2 - 2 * r * rd * cos(phid))
    common_projection = 1.0 / (rd^2 * delta_sq^2)
    v[:] = nt * common_projection * [r_projection, z]
end

"""
Integrates the radiation acceleration integrand on the disc.

# Parameters
- r : radial coordinate [Rg]
- z : height coordinate [Rg]
-
"""
function integrate_radiation_force_integrand(
    radiative_transfer::RadiativeTransfer,
    r,
    z,
    rmin,
    rmax = 1600;
    phi_min = 0.0,
    phi_max = π,
    atol = 0,
    rtol = 1e-4,
    norm = Cubature.INDIVIDUAL,
    maxevals = 0,
)
    f(x, v) = radiation_force_integrand!(radiative_transfer, v, x[1], x[2], r, z)
    return hcubature(
        2,
        f,
        (rmin, phi_min),
        (rmax, phi_max),
        abstol = atol,
        reltol = rtol,
        error_norm = norm,
        maxevals = maxevals,
    )
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
    radiative_transfer::RERadiativeTransfer,
    r,
    z,
    rmax = 1600,
    atol = 0.0,
    rtol = 1e-4,
    norm = Cubature.INDIVIDUAL,
    maxevals = 0,
)
    #println("r : $r, z : $z")
    res, err = integrate_radiation_force_integrand(
        radiative_transfer,
        r,
        z,
        radiative_transfer.radiation.isco,
        rmax,
        atol = atol,
        rtol = rtol,
        norm = norm,
        maxevals = maxevals,
    )
    radiation_constant = compute_radiation_constant(radiative_transfer.radiation)
    force = z * radiation_constant .* res
    #println("force $force")
    return force
end
