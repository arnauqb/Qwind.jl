using Cubature, TimerOutputs
export compute_disc_radiation_field


function compute_tau_uv_integrand(
    radiation::Radiation,
    tau_uv_calculation::TauUVDisk;
    rd,
    phid,
    r,
    z,
)
    return compute_tau_uv(radiation, rd = rd, phid = phid, r = r, z = z)
end

compute_tau_uv_integrand(
    radiation::Radiation,
    tau_uv_calculation::Union{TauUVCenter,NoTauUV};
    rd,
    phid,
    r,
    z,
) = 0.0

compute_tau_uv_integrand(radiation; rd, phid, r, z) = compute_tau_uv_integrand(
    radiation,
    radiation.tau_uv_calculation,
    rd = rd,
    phid = phid,
    r = r,
    z = z,
)

function radiation_force_integrand!(
    v,
    radiation::Radiation;
    rd,
    phid,
    r,
    z,
    vr,
    vz,
    beta,
    gamma,
)
    delta = distance_from_disk(rd, phid, radiation.z_disk, r, 0.0, z)
    r_projection = (r - rd * cos(phid))
    common_projection = 1.0 / (rd^2 * delta^4)
    nt = disk_nt_rel_factors(radiation, rd)
    rel = relativistic_correction(
        radiation.relativistic,
        r = r,
        z = z,
        z0 = radiation.z_disk,
        vr = vr,
        vz = vz,
        beta = beta,
        gamma = gamma,
        r_projection = r_projection,
        delta = delta,
    )
    tau_uv = compute_tau_uv_integrand(radiation, rd = rd, phid = phid, r = r, z = z)
    fuv, mdot = get_fuv_mdot(radiation, rd)
    factors = fuv * mdot * rel * nt * exp(-tau_uv) * common_projection
    v[1] = factors * r_projection
    v[2] = factors * (z - radiation.z_disk)
end

force_prefactors(radiation::Radiation, flag::TauUVCalculationFlag, r, z) =
    3 / (π * radiation.bh.efficiency)
function force_prefactors(radiation::Radiation, ::TauUVCenter, r, z)
    constant = 3 / (π * radiation.bh.efficiency)
    tau_uv = compute_uv_tau(radiation, rd = 0.0, phid = 0.0, r = r, z = z)
    d = sqrt(r^2 + z^2)
    return constant * exp(-tau_uv .* [r / d, z / d])
end
force_prefactors(radiation::Radiation, r, z) =
    force_prefactors(radiation, radiation.tau_uv_calculation, r, z)

"""
Integrates the radiation acceleration integrand on the disc.

# Parameters
- r : radial coordinate [Rg]
- z : height coordinate [Rg]
-
"""
function integrate_radiation_force_integrand(
    radiation;
    r,
    z,
    vr,
    vz,
    rmin,
    rmax = 1600,
    phi_min = 0.0,
    phi_max = π,
    atol = 0.0,
    rtol = 1e-4,
    norm = Cubature.INDIVIDUAL,
    maxevals = 50000,
)
    beta = compute_beta(vr, vz)
    gamma = compute_gamma(beta)
    f(x, v) = radiation_force_integrand!(
        v,
        radiation,
        rd = x[1],
        phid = x[2],
        r = r,
        z = z,
        vr = vr,
        vz = vz,
        beta = beta,
        gamma = gamma,
    )
    integrated, error = hcubature(
        2,
        f,
        (rmin, phi_min),
        (rmax, phi_max),
        abstol = atol,
        reltol = rtol,
        error_norm = norm,
        maxevals = maxevals,
    )
    return integrated
end

function compute_disc_radiation_field_small_heights(radiation::Radiation; r, z, vr, vz)
    Rg = radiation.bh.Rg
    constant = 3 / (2 * radiation.bh.efficiency)
    fuv, mdot = get_fuv_mdot(radiation, r)
    tau_uv = compute_tau_uv(radiation, rd = r, phid = 0.0, r = r, z = z)
    beta = compute_beta(vr, vz)
    rel = ((1 - beta) / (1 + beta))^2
    nt = disk_nt_rel_factors(radiation, r)
    return [0.0, nt * rel * constant * mdot * exp(-tau_uv) * fuv / r^3]
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
    radiation::Radiation;
    r,
    z,
    vr,
    vz,
    rmax = 1600,
    atol = 0,
    rtol = 1e-3,
    norm = Cubature.INDIVIDUAL,
    maxevals = 10000,
    max_z_vertical_flux = 5e-1,
)
    if (z - radiation.z_disk) < max_z_vertical_flux
        if r < radiation.disk_r_in
            return [0.0, 0.0]
        end
        force = compute_disc_radiation_field_small_heights(
            radiation,
            r = r,
            z = z,
            vr = vr,
            vz = vz,
        )
    else
        integration = integrate_radiation_force_integrand(
            radiation,
            r = r,
            z = z,
            vr = vr,
            vz = vz,
            rmin = radiation.disk_r_in,
            rmax = rmax,
            atol = atol,
            rtol = rtol,
            norm = norm,
            maxevals = maxevals,
        )
        prefactors = force_prefactors(radiation, r, z)
        force = (z - radiation.z_disk) .* prefactors .* integration
    end
    return force
end
