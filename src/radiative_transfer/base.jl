export compute_disc_radiation_field

"""
Radiation force integrand for the RE model.
"""
function radiation_force_integrand!(
    radiative_transfer::RadiativeTransfer,
    integration_type::IntegrationFromCenter,
    v,
    rd,
    phid,
    r,
    z,
)
    nt = disk_nt_rel_factors(radiative_transfer.radiation, rd)
    tauuv = compute_uv_tau(
        radiative_transfer,
        rd,
        r,
        z,
    )
    fuv, mdot = get_fuv_mdot(radiative_transfer.radiation, rd)
    r_projection = (r - rd * cos(phid))
    delta_sq = (r^2 + rd^2 + z^2 - 2 * r * rd * cos(phid))
    common_projection = 1.0 / (rd^2 * delta_sq^2)
    v[:] = exp(-tauuv) * fuv * mdot * nt * common_projection * [r_projection, z]
end

function radiation_force_integrand!(
    rt::RadiativeTransfer,
    integration_type::IntegrationFromStreamline,
    v,
    p,
    psi,
    r,
    z,
)
    cosψ = cos(psi)
    rd = sqrt(r^2 + p^2 - 2 * p * r * cosψ)
    if rd <= 6.0 || rd > 1600.0
        v[1] = 0.0
        v[2] = 0.0
        return
    end
    nt = disk_nt_rel_factors(rt.radiation, rd)
    fuv, mdot = get_fuv_mdot(rt.radiation, rd)
    delta2 = p^2 + z^2
    tauuv = compute_uv_tau(rt, rd, r, z)
    #println("rd $rd r $r z $z tau_uv $tauuv")
    tauuv = tauuv / sqrt((r - rd)^2 + z^2) * sqrt(delta2)
    common_part = nt * mdot * fuv * p / delta2^2 / rd^3 * exp(-tauuv)
    v[1] = common_part * p * cosψ
    v[2] = common_part
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
    maxevals = 50000,
    zmax_fromstreamline = 1e-4,
)
    if z < zmax_fromstreamline
        #return [0.0, 0.0], [0.0, 0.0]
        integration_type = IntegrationFromStreamline()
    else
        integration_type = IntegrationFromCenter()
    end
    f(x, v) = radiation_force_integrand!(
        radiative_transfer,
        integration_type,
        v,
        x[1],
        x[2],
        r,
        z,
    )
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
    radiative_transfer::RadiativeTransfer,
    r,
    z;
    rmax = 1600,
    atol = 0,
    rtol = 1e-2,
    norm = Cubature.INDIVIDUAL,
    maxevals = 10000,
)
    #println("r : $r,\t z : $z")
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
