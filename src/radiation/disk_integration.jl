using Cubature, TimerOutputs
export compute_disc_radiation_field


function compute_tau_uv_integrand(
    radiation::Radiation,
    tau_uv_calculation::TauUVDisk;
    r_disk,
    phi_disk,
    r_wind,
    z_wind,
    z_disk,
)
    return compute_tau_uv(
        radiation,
        ri = r_disk,
        phii = phi_disk,
        zi = z_disk,
        rf = r_wind,
        phif = 0.0,
        zf = z_wind,
    )
end

compute_tau_uv_integrand(
    radiation::Radiation,
    tau_uv_calculation::Union{TauUVCenter,NoTauUV};
    r_disk,
    phi_disk,
    r_wind,
    z_wind,
    z_disk,
) = 0.0

compute_tau_uv_integrand(
    radiation,
    tau_uv_calculation;
    r_disk,
    phi_disk,
    r_wind,
    z_wind,
    z_disk,
) = compute_tau_uv_integrand(
    radiation,
    tau_uv_calculation,
    r_disk = r_disk,
    phi_disk = phi_disk,
    r_wind = r_wind,
    z_wind = z_wind,
    z_disk = z_disk,
)

function radiation_force_integrand!(
    v,
    radiation::Radiation,
    parameters::Parameters;
    r_disk,
    phi_disk,
    r_wind,
    z_wind,
    vr_wind,
    vz_wind,
    beta,
    gamma,
)
    delta = compute_distance_cylindrical(
        r_disk,
        phi_disk,
        parameters.z_disk,
        r_wind,
        0.0,
        z_wind,
    )
    r_projection = (r_wind - r_disk * cos(phi_disk))
    common_projection = 1.0 / (r_disk^2 * delta^4)
    nt = disk_nt_rel_factors(radiation, r_disk)
    rel = relativistic_correction(
        parameters.relativistic_flag,
        r = r_wind,
        z = z_wind,
        z0 = parameters.z_disk,
        vr = vr_wind,
        vz = vz_wind,
        beta = beta,
        gamma = gamma,
        r_projection = r_projection,
        delta = delta,
    )
    tau_uv = compute_tau_uv_integrand(
        radiation,
        parameters.tau_uv_calculation_flag,
        r_disk = r_disk,
        phi_disk = phi_disk,
        r_wind = r_wind,
        z_wind = z_wind,
        z_disk = parameters.z_disk,
    )
    fuv, mdot = get_fuv_mdot(radiation, r_disk)
    factors = fuv * mdot * rel * nt * exp(-tau_uv) * common_projection
    v[1] = factors * r_projection
    v[2] = factors * (z_wind - parameters.z_disk)
end

force_prefactors(radiation::Radiation, flag::TauUVCalculationFlag, r_wind, z_wind, z_disk) =
    3 / (π * radiation.bh.efficiency)

function force_prefactors(radiation::Radiation, ::TauUVCenter, r_wind, z_wind, z_disk)
    constant = 3 / (π * radiation.bh.efficiency)
    tau_uv = compute_tau_uv(
        radiation,
        ri = 0.0,
        phii = 0.0,
        zi = z_disk,
        rf = r_wind,
        phif = 0.0,
        zf = z_wind,
    )
    d = sqrt(r_wind^2 + z_wind^2)
    return constant * exp.(-tau_uv)# .* [r / d, z / d])
end

"""
Integrates the radiation acceleration integrand on the disc.

# Parameters
- r : radial coor_diskinate [Rg]
- z : height coor_diskinate [Rg]
-
"""
function integrate_radiation_force_integrand(
    radiation,
    parameters;
    r_wind,
    z_wind,
    vr_wind,
    vz_wind,
    norm = Cubature.INDIVIDUAL,
    maxevals = 50000,
    disk_r_in,
    disk_r_out,
)
    beta = compute_beta(vr_wind, vz_wind)
    gamma = compute_gamma(beta)
    f(x, v) = radiation_force_integrand!(
        v,
        radiation,
        parameters,
        r_disk = x[1],
        phi_disk = x[2],
        r_wind = r_wind,
        z_wind = z_wind,
        vr_wind = vr_wind,
        vz_wind = vz_wind,
        beta = beta,
        gamma = gamma,
    )
    integrated, error = hcubature(
        2,
        f,
        (disk_r_in, 0.0),
        (disk_r_out, π),
        abstol = parameters.disk_integral_atol,
        reltol = parameters.disk_integral_rtol,
        error_norm = norm,
        maxevals = maxevals,
    )
    return integrated
end

function compute_disc_radiation_field_small_heights(
    radiation::Radiation,
    parameters::Parameters;
    r_wind,
    z_wind,
    vr_wind,
    vz_wind,
)
    constant = 3 / (2 * radiation.bh.efficiency)
    fuv, mdot = get_fuv_mdot(radiation, r_wind)
    tau_uv = compute_tau_uv_integrand(
        radiation,
        parameters.tau_uv_calculation_flag,
        r_disk = r_wind,
        phi_disk = 0.0,
        r_wind = r_wind,
        z_wind = z_wind,
        z_disk = parameters.z_disk,
    )
    beta = compute_beta(vr_wind, vz_wind)
    rel = ((1 - beta) / (1 + beta))^2
    nt = disk_nt_rel_factors(radiation, r_wind)
    return [0.0, nt * rel * constant * mdot * exp(-tau_uv) * fuv / r_wind^3]
end

"""
Computes the disc radiation field at the point (r,z) by performing
an integral over the disc.

# Parameters
- r: radial coor_diskinate [Rg]
- z: height coor_diskinate [Rg]
- radiative_efficiency: accretion radiative efficiency
"""
function compute_disc_radiation_field(
    radiation::Radiation,
    parameters::Parameters;
    r_wind,
    z_wind,
    vr_wind,
    vz_wind,
    norm = Cubature.INDIVIDUAL,
    maxevals = 10000,
    max_z_vertical_flux = 5e-1,
)
    if parameters.disk_r_in == "corona_radius"
        disk_r_in = radiation.qsosed_model.corona.radius
    elseif parameters.disk_r_in == "warm_radius"
        disk_r_in = radiation.qsosed_model.warm.radius
    else
        disk_r_in = parameters.disk_r_in
    end
    if (z_wind - parameters.z_disk) < max_z_vertical_flux
        if r_wind < disk_r_in
            return [0.0, 0.0]
        end
        force = compute_disc_radiation_field_small_heights(
            radiation,
            parameters,
            r_wind = r_wind,
            z_wind = z_wind,
            vr_wind = vr_wind,
            vz_wind = vz_wind,
        )
    else
        integration = integrate_radiation_force_integrand(
            radiation,
            parameters,
            r_wind = r_wind,
            z_wind = z_wind,
            vr_wind = vr_wind,
            vz_wind = vz_wind,
            disk_r_in = disk_r_in,
            disk_r_out = parameters.disk_r_out,
        )
        prefactors = force_prefactors(
            radiation,
            parameters.tau_uv_calculation_flag,
            r_wind,
            z_wind,
            parameters.z_disk,
        )
        force = (z_wind - parameters.z_disk) .* prefactors .* integration
    end
    return force
end
