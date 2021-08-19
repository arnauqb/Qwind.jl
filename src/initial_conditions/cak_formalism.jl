using Peaks, Roots
export calculate_wind_mdots

get_B0(r) = 1.0 / r^2

function gamma0(radiation, r)
    ret = compute_disc_radiation_field_small_heights(
        radiation,
        r = r,
        z = 1.0,
        vr = 0.0,
        vz = 0.0,
    )[2]
    # remove attenuation and fuv
    fuv, mdot = get_fuv_mdot(radiation, r)
    tau_uv = compute_tau_uv(radiation, rd = r, phid = 0.0, r = r, z = 1.0)
    ret = ret / fuv / exp(-tau_uv)
    return ret
end

function f(radiation::Radiation, z; r, alpha = 0.6, zmax = 5e-1)
    cc = 1 / (alpha^alpha * (1 - alpha)^(1 - alpha))
    frad = compute_disc_radiation_field(
        radiation,
        r = r,
        z = z,
        vr = 0.0,
        vz = 0.0,
        max_z_vertical_flux = zmax,
        rtol = 1e-4,
        maxevals = 100000,
    )[2]
    # Restore
    radiation.fuv_grid .= fuv_copy
    radiation.wi.density_grid.grid .= dgrid_copy
    B0 = get_B0(r)
    return -(grav + fr) / B0
end

function nozzle_function(radiation::Radiation, z; r, alpha = 0.6, zmax = 5e-1)
    c = alpha * (1 - alpha)^((1 - alpha) / alpha)
    if g(radiation, z, r = r, zmax = zmax) <= 0
        return Inf
    end
    return c * f(radiation, z, r = r, alpha = alpha, zmax = zmax)^(1 / alpha) /
           (g(radiation, z, r = r, zmax = zmax))^((1 - alpha) / alpha)
end

function find_nozzle_function_minimum(
    radiation::Radiation,
    r;
    alpha = 0.6,
    zmax = 5e-1,
    n_z = 150,
)
    z_range = 10 .^ range(-4, 3, length = n_z)
    n_range =
        nozzle_function.(
            Ref(radiation),
            z_range,
            r = r,
            alpha = alpha,
            zmax = zmax,
        )
    # check if it's monotonic
    mon_increasing = reduce(*, aaccumulate(max, n_range) .== n_range)
    mon_decreasing = reduce(*, aaccumulate(min, n_range) .== n_range)
    if mon_increasing || mon_decreasing
        return Inf, NaN
    end
    mask = n_range .!= Inf
    n_range = n_range[mask]
    z_range = z_range[mask]
    minima_arg, minima_values = findminima(n_range)
    if length(minima_values) == 0
        return Inf, NaN
    end
    minn = argmin(minima_values)
    minnarg = minima_arg[minn]
    return z_range[minnarg], minima_values[minn]
end

function CAK_Σ(radiation::Radiation, r; K = 0.3, alpha = 0.6, mu_e = 1.17)
    cc = alpha * (1 - alpha)^((1 - alpha) / alpha)
    vth = compute_thermal_velocity(25e3) * C
    B0 = get_B0(r) * C^2 / radiation.bh.Rg
    γ0 = gamma0(radiation, r)
    constant = K * (1 / (SIGMA_E / mu_e * vth))^alpha * C^2 / radiation.bh.Rg
    return cc * constant * γ0^(1 / alpha) / B0^((1 - alpha) / alpha)
end

function get_initial_density(radiation::Radiation; r, mdot, K = 0.3, alpha = 0.6, mu = 0.61)
    sigmadot = mdot * CAK_Σ(radiation, r, K = K, alpha = alpha)
    return sigmadot / (compute_thermal_velocity(disk_temperature(radiation.bh, r)) * C) /
           (mu * M_P)
end

function calculate_wind_mdots(radiation::Radiation)
    rr = 10 .^ range(log10(15.0), log10(1500), length = 50)
    mdots = []
    zcs = []
    f(r) = find_nozzle_function_minimum(radiation, r, alpha = 0.6, zmax = 1e-2)
    results = pmap(f, rr)
    zcs = [res[1] for res in results]
    mdots = [res[2] for res in results]
    rr = rr[.!isnan.(mdots)]
    mdots = mdots[.!isnan.(mdots)]
    zcs = zcs[zcs .!= Inf]
    return rr, mdots, zcs
end
