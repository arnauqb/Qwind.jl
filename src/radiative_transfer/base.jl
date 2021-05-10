export compute_disc_radiation_field

function radiation_force_integrand!(
    flux_correction::NoRelativistic,
    radiative_transfer::RadiativeTransfer,
    radiation::Radiation,
    density_grid::InterpolationGrid,
    Rg,
    v,
    rd,
    phid,
    r,
    z,
    vr,
    vz,
    beta,
    gamma,
)
    delta = sqrt(r^2 + rd^2 + z^2 - 2 * r * rd * cos(phid))
    nt = disk_nt_rel_factors(radiation, rd)
    tauuv = compute_uv_tau(density_grid, density_grid.iterator, rd, 0.0, r, z, Rg)
    # deproject tauuv
    tauuv = tauuv * delta / d_euclidean(rd, r, 0.0, z)
    fuv, mdot = get_fuv_mdot(radiation, rd)
    r_projection = (r - rd * cos(phid))
    common_projection = 1.0 / (rd^2 * delta^4)
    v[:] = exp(-tauuv) * fuv * mdot * nt * common_projection * [r_projection, z]
end

# Relativistic version
function radiation_force_integrand!(
    flux_correction::Relativistic,
    radiative_transfer::RadiativeTransfer,
    radiation::Radiation,
    density_grid::InterpolationGrid,
    Rg,
    v,
    rd,
    phid,
    r,
    z,
    vr,
    vz,
    beta,
    gamma,
)
    delta = sqrt(r^2 + rd^2 + z^2 - 2 * r * rd * cos(phid))
    tauuv = compute_uv_tau(density_grid, density_grid.iterator, rd, 0.0, r, z, Rg)
    # deproject tauuv
    tauuv = tauuv * delta / d_euclidean(rd, r, 0.0, z)
    fuv, mdot = get_fuv_mdot(radiation, rd)
    nt = disk_nt_rel_factors(r, radiation.spin, radiation.isco)
    r_projection = (r - rd * cos(phid))
    # relativistic correction to the flux
    cosθ = (r_projection * vr + z * vz) / (delta * beta)
    flux_correction = 1.0 / (gamma * (1 + beta * cosθ))^4
    # common geometric term for r and z
    common_projection = 1.0 / (rd^2 * delta^4)
    v[:] =
        exp(-tauuv) *
        flux_correction *
        fuv *
        mdot *
        nt *
        common_projection *
        [r_projection, z]
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
    vr,
    vz,
    beta,
    gamma,
    rmin,
    rmax = 1600;
    phi_min = 0.0,
    phi_max = π,
    atol = 0.0,
    rtol = 1e-4,
    norm = Cubature.INDIVIDUAL,
    maxevals = 50000,
)
    f(x, v) = radiation_force_integrand!(
        radiative_transfer.radiation.flux_correction,
        radiative_transfer,
        radiative_transfer.radiation,
        radiative_transfer.interpolator.density_grid,
        radiative_transfer.radiation.Rg,
        v,
        x[1],
        x[2],
        r,
        z,
        vr,
        vz,
        beta,
        gamma,
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


function compute_disc_radiation_field_vertical(
    radiation::Radiation,
    grid::InterpolationGrid,
    r,
    z,
)
    Rg = radiation.Rg
    constant = 3 / (2 * radiation.efficiency)
    fuv, mdot = get_fuv_mdot(radiation, r)
    tauuv = compute_uv_tau(grid, r, r, z, Rg)
    nt = disk_nt_rel_factors(r, radiation.spin, radiation.isco)
    return [0.0, nt * constant * mdot * exp(-tauuv) * fuv / r^3]
end

function compute_disc_radiation_field_vertical(
    flux_correction::Relativistic,
    radiation::Radiation,
    grid::InterpolationGrid,
    r,
    z,
    beta,
)
    flux_correction = ((1 - beta) / (1 + beta))^2
    return flux_correction * compute_disc_radiation_field_vertical(radiation, grid, r, z)
end

compute_disc_radiation_field_vertical(
    flux_correction::NoRelativistic,
    radiation::Radiation,
    grid::InterpolationGrid,
    r,
    z,
    beta,
) = compute_disc_radiation_field_vertical(radiation, grid, r, z)

compute_disc_radiation_field_vertical(radiative_transfer::RadiativeTransfer, r, z, beta) =
    compute_disc_radiation_field_vertical(
        radiative_transfer.radiation.flux_correction,
        radiative_transfer.radiation,
        radiative_transfer.interpolator.density_grid,
        r,
        z,
        beta,
    )

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
    z,
    vr,
    vz;
    rmax = 1600,
    atol = 0,
    rtol = 1e-3,
    norm = Cubature.INDIVIDUAL,
    maxevals = 10000,
    max_z_vertical_flux = 5e-1,#1e4 #5e-1,
)
    # sometimes the solver may try unphysical values of vr and vz, so we take the max
    # Similarly for the lower bound as beta = 0 leads to numerical singularity
    beta = max(min(1.0, sqrt(vr^2 + vz^2)), 1e-4)
    gamma = 1.0 / sqrt(1 - beta^2)
    if z < max_z_vertical_flux
        if r < radiative_transfer.radiation.disk_r_in
            return [0.0, 0.0]
        end
        force = compute_disc_radiation_field_vertical(radiative_transfer, r, z, beta)
    else
        res, err = integrate_radiation_force_integrand(
            radiative_transfer,
            r,
            z,
            vr,
            vz,
            beta,
            gamma,
            radiative_transfer.radiation.disk_r_in,
            rmax,
            atol = atol,
            rtol = rtol,
            norm = norm,
            maxevals = maxevals,
        )
        radiation_constant = compute_radiation_constant(radiative_transfer.radiation)
        force = z * radiation_constant .* res
    end
    return force
end
