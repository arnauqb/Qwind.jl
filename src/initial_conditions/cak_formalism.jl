using Peaks, Roots
export calculate_wind_mdots

get_B0(r) = 1.0 / r^2

function gamma0(density_grid, iterator, parameters, radiation, r)
    parameters_no_tau = change_parameter(parameters, :tau_uv_calculation_flag, NoTauUV())
    ret = compute_disc_radiation_field_small_heights(
        density_grid,
        iterator,
        parameters_no_tau,
        radiation,
        r_wind = r,
        z_wind = 1.0,
        vr_wind = 0.0,
        vz_wind = 0.0,
    )[2]
    fuv, mdot = get_fuv_mdot(radiation, r)
    # remove fuv
    ret = ret / fuv
    return ret
end

function f(density_grid, iterator, parameters, radiation, z; r, alpha = 0.6, zmax = 1e-1)
    cc = 1 / (alpha^alpha * (1 - alpha)^(1 - alpha))
    parameters_no_tau = change_parameter(parameters, :tau_uv_calculation_flag, NoTauUV())
    frad = compute_disc_radiation_field(
        density_grid,
        iterator,
        parameters_no_tau,
        radiation,
        r_wind = r,
        z_wind = z,
        vr_wind = 0.0,
        vz_wind = 0.0,
        max_z_vertical_flux = zmax,
        maxevals = 100000,
    )[2]
    f0 = gamma0(density_grid, iterator, parameters, radiation, r)
    return cc * frad / f0
end

function g(density_grid, iterator, parameters, radiation, z; r, zmax = 1e-1)
    if parameters.ic_z_disk == "auto"
        z_disk = disk_height(radiation.bh, r)
    else
        z_disk = parameters.ic_z_disk
    end
    grav = compute_gravitational_acceleration(r, z + z_disk)[2]
    fuv_copy = copy(radiation.fuv_grid)
    parameters_no_tau = change_parameter(parameters, :tau_uv_calculation_flag, NoTauUV())
    radiation.fuv_grid .= 1.0
    fr = compute_disc_radiation_field(
        density_grid,
        iterator,
        parameters_no_tau,
        radiation,
        r_wind = r,
        z_wind = z,
        vr_wind = 0.0,
        vz_wind = 0.0,
        max_z_vertical_flux = zmax,
        maxevals = 100000,
    )[2]
    B0 = get_B0(r)
    T = disk_temperature(radiation.bh, r)
    vth = compute_thermal_velocity(T)
    s = vth^2 / (2 * B0 * r)
    a = 1 + (z / r)^2
    gas_pressure_term = 4 * s * z / r / a
    radiation.fuv_grid .= fuv_copy
    return -(grav + fr) / B0 + gas_pressure_term
end

function nozzle_function(
    density_grid,
    iterator,
    parameters,
    radiation,
    z;
    r,
    alpha = 0.6,
    zmax = 1e-1,
    include_a = true,
)
    c = alpha * (1 - alpha)^((1 - alpha) / alpha)
    if include_a
        a = 1 + (z / r)^2
    else
        a = 1.0
    end
    if g(density_grid, iterator, parameters, radiation, z, r = r, zmax = zmax) <= 0
        return Inf
    end
    return c *
           a *
           f(
               density_grid,
               iterator,
               parameters,
               radiation,
               z,
               r = r,
               alpha = alpha,
               zmax = zmax,
           )^(1 / alpha) / (
        g(
            density_grid,
            iterator,
            parameters,
            radiation,
            z,
            r = r,
            zmax = zmax,
        )^((1 - alpha) / alpha)
    )
end

function find_nozzle_function_minimum(
    density_grid,
    iterator,
    parameters,
    radiation,
    r;
    alpha = 0.6,
    zmax = 1e-1,
    nz = 250,
    include_a = true,
)
    if zmax == Inf
        z0 = -2
    else
        z0 = log10(zmax)
    end
    z0 = max(-2, log10(parameters.z_disk))
    z_range = 10 .^ range(z0, 3, length = nz)
    n_range =
        nozzle_function.(
            Ref(density_grid),
            Ref(iterator),
            Ref(parameters),
            Ref(radiation),
            z_range,
            r = r,
            alpha = alpha,
            zmax = zmax,
            include_a = include_a,
        )
    # check if it's monotonic
    mon_increasing = reduce(*, accumulate(max, n_range) .== n_range)
    mon_decreasing = reduce(*, accumulate(min, n_range) .== n_range)
    if mon_increasing || mon_decreasing
        return Inf, NaN
    end
    minv, mina = findmin(n_range)
    return z_range[mina], minv
end

function CAK_Σ(density_grid, iterator, parameters, radiation, r; K = nothing)
    if K == nothing
        K = parameters.ic_K
    end
    alpha = parameters.ic_alpha
    cc = alpha * (1 - alpha)^((1 - alpha) / alpha)
    vth = compute_thermal_velocity(25e3) * C
    B0 = get_B0(r) * C^2 / radiation.bh.Rg
    constant = K * (1 / (SIGMA_E * vth))^alpha * C^2 / radiation.bh.Rg
    γ0 = gamma0(density_grid, iterator, parameters, radiation, r) * constant
    return cc * γ0^(1 / alpha) / B0^((1 - alpha) / alpha)
end

function get_initial_density(radiation::Radiation, wind, parameters; r, mdot, zc)
    if parameters.ic_K == "auto"
        tau_x = compute_tau_xray(
            wind.density_grid,
            wind.grid_iterator,
            radiation,
            parameters,
            ri = 0.0,
            phii = 0.0,
            zi = parameters.z_xray,
            rf = r,
            phif = 0.0,
            zf = zc,
        )
        n = interpolate_density(wind.density_grid, r, zc)
        xi = compute_ionization_parameter(
            radiation,
            wind,
            parameters,
            r = r,
            z = zc,
            number_density = n,
            tau_x = tau_x,
        )
        K = compute_force_multiplier_k(xi, FMNoInterp())
    else
        K = parameters.ic_K
    end
    sigmadot =
        mdot * CAK_Σ(wind.density_grid, wind.grid_iterator, parameters, radiation, r, K = K)
    density =
        sigmadot / (
            compute_thermal_velocity(
                disk_temperature(radiation.bh, r),
                parameters.mu_nucleon,
            ) * C
        )
    number_density = density / M_P / parameters.mu_nucleon
    return number_density
end

function calculate_wind_mdots(
    density_grid,
    iterator,
    parameters,
    radiation;
    rmin = 6.1,
    rmax = 1500.0,
    nr = 100,
    nz = 100,
)
    rr = 10 .^ range(log10(rmin), log10(rmax), length = nr)
    mdots = []
    zcs = []
    f(r) = find_nozzle_function_minimum(
        density_grid,
        iterator,
        parameters,
        radiation,
        r,
        alpha = 0.6,
        zmax = 1e-1,
        nz = nz,
    )
    results = pmap(f, rr)
    zcs = [res[1] for res in results]
    mdots = [res[2] for res in results]
    rr = rr[.!isnan.(mdots)]
    zcs = zcs[.!isnan.(mdots)]
    mdots = mdots[.!isnan.(mdots)]
    return rr, mdots, zcs
end
