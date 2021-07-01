using Peaks, Roots

get_B0(bh, r) = 1.0 / r^2

function f(rt::RadiativeTransfer, bh::BlackHole, z; r, alpha = 0.6, zmax = 5e-1)
    cc = 1 / (alpha^alpha * (1 - alpha)^(1 - alpha))
    frad = compute_disc_radiation_field(
        rt,
        r,
        z,
        0.0,
        0.0,
        max_z_vertical_flux = zmax,
        rtol = 1e-4,
        maxevals = 100000,
    )[2]
    frad0 = compute_disc_radiation_field(
        rt,
        r,
        z,
        0.0,
        0.0,
        max_z_vertical_flux = Inf,
        no_uv_fraction = true,
    )[2]
    return cc * frad / frad0
end

function g(rt::RadiativeTransfer, bh::BlackHole, z; r, zmax = 5e-1)
    grav = compute_gravitational_acceleration(r, z + disk_height(bh, r))[2]
    fr = compute_disc_radiation_field(
        rt,
        r,
        z,
        0.0,
        0.0,
        max_z_vertical_flux = zmax,
        rtol = 1e-4,
        maxevals = 100000,
        no_uv_fraction = true,
    )[2]
    B0 = get_B0(bh, r)
    return -(grav + fr) / B0
end

function nozzle_function(
    rt::RadiativeTransfer,
    bh::BlackHole,
    z;
    r,
    alpha = 0.6,
    zmax = 5e-1,
)
    c = alpha * (1 - alpha)^((1 - alpha) / alpha)
    if g(rt, bh, z, r = r, zmax = zmax) <= 0
        return Inf
    end
    return c * f(rt, bh, z, r = r, alpha = alpha, zmax = zmax)^(1 / alpha) /
           (g(rt, bh, z, r = r, zmax = zmax))^((1 - alpha) / alpha)
end

function find_nozzle_function_minimum(
    rt::RadiativeTransfer,
    bh::BlackHole,
    r;
    alpha = 0.6,
    zmax = 5e-1,
)
    z_range = 10 .^ range(-6, 6, length = 500)
    n_range = nozzle_function.(Ref(rt), Ref(bh), z_range, r = r, alpha = alpha, zmax = zmax)
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

function beta(bh::BlackHole, vz; r)
    cs = compute_thermal_velocity(disk_temperature(bh, r))
    return 1 - cs^2 / vz^2
end

function CAK_Σ(rt::RadiativeTransfer, bh::BlackHole, r; K = 0.3, alpha = 0.6)
    cc = alpha * (1 - alpha)^((1 - alpha) / alpha)
    vth = compute_thermal_velocity(25e3) * C
    B0 = get_B0(bh, r) * C^2 / bh.Rg
    gamma0 =
        compute_disc_radiation_field(
            rt,
            r,
            r,
            0.0,
            0.0,
            max_z_vertical_flux = Inf,
            no_uv_fraction = true,
        )[2] *
        K *
        (1 / (SIGMA_E * vth))^alpha *
        C^2 / bh.Rg
    return cc * gamma0^(1 / alpha) / B0^((1 - alpha) / alpha)
end

function get_wp0(rt::RadiativeTransfer, bh::BlackHole, z, w; r, mdot, alpha = 0.6)
#    wp0 = (g(rt, bh, z, r = r) / f(rt, bh, z, r = r, alpha = alpha))^(1 / alpha) * mdot
    wp0 = mdot * (1 / f(rt, bh, z, r=r, alpha=alpha))^(1/(alpha-1))
    return wp0
end

function Feq(rt::RadiativeTransfer, bh::BlackHole, z, w, wp, mdot; r, alpha = 0.6)
    gv = g(rt, bh, z, r = r)
    fv = f(rt, bh, z, r = r, alpha = alpha)
    vz = sqrt(2 * r * get_B0(bh, r) * w)
    cs = compute_thermal_velocity(disk_temperature(bh, r))
    #ret = wp * (1 - cs^2 / (2w)) - gv - fv * (1 / mdot)^alpha * abs(wp)^alpha
    ret = wp * (1 - cs^2 / (vz^2)) + gv - fv * (1 / mdot)^alpha * abs(wp)^alpha
    return ret
end

function ic_residual!(
    rt::RadiativeTransfer,
    bh::BlackHole,
    out,
    du,
    u,
    p,
    t;
    r,
    mdot,
    alpha = 0.6,
)
    z = t
    w = u[1]
    wp = du[1]
    out[1] = Feq(rt, bh, z, w, wp, mdot, r = r, alpha = alpha)
end

function initialize_ic_integrator(
    rt::RadiativeTransfer,
    bh::BlackHole,
    mdot,
    r;
    alpha = 0.6,
)
    cs = compute_thermal_velocity(disk_temperature(bh, r))
    z0 = disk_height(bh, r) / 2
    u0 = [cs^2 / (2 * get_B0(bh, r) * r)]
    du0 = [get_wp0(rt, bh, z0, u0[1], r = r, mdot = mdot, alpha = alpha)]
    z_range = (z0, 30)
    dae_residual!(out, du, u, p, t) =
        ic_residual!(rt, bh, out, du, u, p, t, r = r, mdot = mdot, alpha = alpha)
    dae_problem = DAEProblem(dae_residual!, du0, u0, z_range, differential_vars = [true])
    integrator = init(dae_problem, IDA(init_all = false))
    integrator.opts.abstol = 1e-8
    integrator.opts.reltol = 1e-3
    return integrator
end

function find_critical_point(rt::RadiativeTransfer, bh::BlackHole, r, mdot; alpha = 0.6)
    integrator = initialize_ic_integrator(rt, bh, mdot, r, alpha = alpha)
    solve!(integrator)
    z_range =
        10 .^
        range(log10(integrator.sol.t[2]), log10(integrator.sol.t[end - 1]), length = 50)
    n_range = nozzle_function.(Ref(rt), Ref(bh), z_range, r = r, alpha = alpha)
    solutions = integrator.sol(z_range)
    u_range = reduce(vcat, solutions.u)
    beta_range = beta.(Ref(bh), u_range, r = r)
    minima_arg, minima_values = findminima(n_range)
    return z_range[minima_arg]
end

function find_critical_point_mdot_kernel(
    rt::RadiativeTransfer,
    bh::BlackHole,
    r,
    mdot;
    alpha = 0.6,
)
    integrator = initialize_ic_integrator(rt, bh, mdot, r, alpha = alpha)
    solve!(integrator)
    z_range =
        10 .^
        range(log10(integrator.sol.t[2]), log10(integrator.sol.t[end - 1]), length = 200)
    zh = Qwind.disk_height(bh, r)
    n_range = nozzle_function.(Ref(rt), Ref(bh), z_range, r = r, alpha = alpha)
    solutions = integrator.sol(z_range)
    u_range = reduce(vcat, solutions.u)
    beta_range = Qwind.beta.(Ref(bh), u_range, r = r)
    minima_arg, minima_values = findminima(n_range)
    if length(minima_arg) == 0
        return Inf, Inf
    end
    return beta_range[minima_arg] * mdot - n_range[minima_arg], z_range[minima_arg]
end

function find_critical_point_mdot(rt::RadiativeTransfer, bh::BlackHole, r; alpha = 0.6)
    f_kernel = mdot -> find_critical_point_mdot_kernel(rt, bh, r, mdot, alpha = alpha)[1]
    zero = find_zero(f_kernel, (1e-5, 1e3), Bisection(), rtol = 1e-3, atol = 0)
    zc = find_critical_point_mdot_kernel(rt, bh, r, zero, alpha = alpha)[2]
    return zero, zc
    #zeros = find_zeros(f_kernel, 1e-8, 1e-2, rtol = 1e-3, atol = 0)
    #println(zeros)
    #return zeros[0]
end

function get_initial_density(rt, bh::BlackHole, r, mdot; K = 0.3, alpha = 0.6, mu = 0.5)
    sigmadot = mdot * CAK_Σ(rt, bh, r, K = K, alpha = alpha)
    return sigmadot / (compute_thermal_velocity(disk_temperature(bh, r)) * C) / (mu * M_P)
end

function find_initial_denisty(model, r)
    zc, mdc = find_nozzle_function_minimum(
        model.rt,
        model.bh,
        r,
        alpha = model.ic.alpha,
        zmax = 5e-1,
    )
    n = get_initial_density(
        model.rt,
        model.bh,
        r,
        mdc,
        K = model.ic.K,
        alpha = model.ic.alpha,
        mu = model.ic.mu,
    )
    return n
end
