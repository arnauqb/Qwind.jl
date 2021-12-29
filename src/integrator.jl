using DiffEqCallbacks, Sundials, Printf
using Statistics: std, mean
using Roots: find_zero, Bisection
using Distributed
using QuadGK
import Qwind.compute_density
import Base.push!
export initialize_integrator,
    run_integrator!,
    run_integrators!,
    run_integrators_parallel!,
    get_initial_radii_and_linewidths,
    compute_density,
    compute_vt,
    compute_d,
    failed,
    IntegratorParameters,
    escaped

Base.length(integrator::Sundials.IDAIntegrator) = length(integrator.p.data[:r])

function make_save_data(trajectory_id = -1)
    ret = Dict{Symbol,Any}(:trajectory_id => trajectory_id)
    keys = [:r, :z, :vr, :vphi, :vz, :n, :fm, :xi, :dvdr, :xi, :taueff, :taux, :tauuv]
    for key in keys
        ret[key] = Float64[]
    end
    return ret
end

struct IntegratorParameters
    id::Int
    r0::Float64
    z0::Float64
    v0::Float64
    n0::Float64
    l0::Float64
    lwnorm::Float64
    r_min::Float64
    r_max::Float64
    z_min::Float64
    z_max::Float64
    data::Any
    finished::Vector{Bool}
    disk_integral_atol::Float64
    disk_integral_rtol::Float64
end

function out_of_grid(p::IntegratorParameters, r, z)
    if p.r_min <= r <= p.r_max && p.z_min <= z <= p.z_max
        return false
    else
        return true
    end
end

function compute_density(r, z, vr, vz, parameters::IntegratorParameters)
    return compute_density(
        r,
        z,
        vr,
        vz,
        parameters.r0,
        parameters.z0,
        parameters.v0,
        parameters.n0,
    )
end

function initialize_integrator(
    rad,
    wind,
    params,
    ic,
    r0,
    linewidth;
    tmax = 1e8,
    trajectory_id = -1,
)
    l0 = getl0(ic, r0)
    z0 = getz0(ic, r0)
    n0 = getn0(ic, radiation = rad, wind = wind, parameters = params, r = r0)
    v0 = getv0(rad.bh, ic, r0, mu_nucleon = params.mu_nucleon)
    lwnorm = linewidth / sqrt(r0^2 + z0^2)
    termination_callback = DiscreteCallback(
        termination_condition,
        integrator -> affect!(integrator),
        save_positions = (false, false),
    )
    data = make_save_data(trajectory_id)
    integ_params = IntegratorParameters(
        trajectory_id,
        r0,
        z0,
        v0,
        n0,
        l0,
        lwnorm,
        params.integrator_r_min,
        params.integrator_r_max,
        params.integrator_z_min,
        params.integrator_z_max,
        data,
        [false],
        params.disk_integral_atol,
        params.disk_integral_rtol,
    )
    saved_data = SavedValues(Float64, Float64)
    saving_callback = SavingCallback(
        (u, t, integrator) -> save(u, t, integrator, rad, wind, params, trajectory_id),
        saved_data,
    )
    callback_set = CallbackSet(termination_callback, saving_callback)
    a₀ = compute_initial_acceleration(rad, wind, params, r0, z0, 0, v0, integ_params)
    du₀ = [0.0, v0, a₀[1], a₀[2]]
    u₀ = [r0, z0, 0.0, v0]
    tspan = (0.0, tmax)
    dae_problem =
        create_dae_problem(rad, wind, params, residual!, du₀, u₀, tspan, integ_params)
    integrator = init(dae_problem, IDA(init_all = false), callback = callback_set)
    integrator.opts.abstol = params.integrator_atol
    integrator.opts.reltol = params.integrator_rtol
    return integrator
end

initialize_integrator(model, r0, linewidth; tmax = 1e8, trajectory_id = -1) =
    initialize_integrator(
        model.rad,
        model.wind,
        model.parameters,
        model.ic,
        r0,
        linewidth,
        tmax = tmax,
        trajectory_id = trajectory_id,
    )

function get_initial_radii_and_linewidths(
    initial_conditions::InitialConditions,
    xray_luminosity,
    Rg,
)
    rin = getrin(initial_conditions)
    rfi = getrfi(initial_conditions)
    nlines = getnlines(initial_conditions)
    if initial_conditions.trajs_spacing == "log"
        lines_widths = diff(10 .^ range(log10(rin), log10(rfi), length = nlines + 1))
        lines_range = [rin + lines_widths[1] / 2.0]
        for i = 2:nlines
            r0 = lines_range[i - 1] + lines_widths[i - 1] / 2 + lines_widths[i] / 2
            push!(lines_range, r0)
        end
    else
        lines_widths = range(rin, rfi, length = nlines)
        dr = (rfi - rin) / nlines
        lines_range = [rin + (i + 0.5) * dr for i = 0:(nlines - 1)]
        lines_widths = diff([lines_range; rfi + dr / 2])
    end
    return lines_range, lines_widths
end

function get_initial_radii_and_linewidths(model)
    nlines = getnlines(model.ic)
    if nlines == "auto"
        lines_range, lines_widths = compute_lines_range(model)
    else
        lines_range, lines_widths = get_initial_radii_and_linewidths(
            model.ic,
            model.rad.xray_luminosity,
            model.bh.Rg,
        )
    end
    return lines_range, lines_widths
end

function create_and_run_integrator(model; r0, linewidth, trajectory_id)
    integrator = initialize_integrator(model, r0, linewidth, trajectory_id = trajectory_id)
    solve!(integrator)
    return integrator
end

run_integrator!(integrator::Sundials.IDAIntegrator) = solve!(integrator)
function run_integrators!(integrators::Vector)
    for (i, integrator) in enumerate(integrators)
        @info "Running integrator $(@sprintf "%03d" i) of $(length(integrators))"
        flush(stdout)
        run_integrator!(integrator)
    end
end

compute_vt(integrator::Sundials.IDAIntegrator) = sqrt(integrator.u[3]^2 + integrator.u[4]^2)

compute_d(integrator::Sundials.IDAIntegrator) = sqrt(integrator.u[1]^2 + integrator.u[2]^2)

out_of_grid(integrator::Sundials.IDAIntegrator) =
    out_of_grid(integrator.p, integrator.u[1], integrator.u[2])

function escaped(integrator::Sundials.IDAIntegrator)
    compute_vt(integrator) > compute_escape_velocity(compute_d(integrator))
end

function failed(integrator::Sundials.IDAIntegrator, r, z)
    intersects = self_intersects(integrator, r, z)
    too_long = length(integrator.p.data[:r]) > 500
    integrator.u[1] < 0.0 ||
        integrator.u[2] < max(integrator.p.z_min, integrator.p.z0) ||
        intersects ||
        too_long
end

compute_density(integrator::Sundials.IDAIntegrator) = compute_density(
    integrator.u[1],
    integrator.u[2],
    integrator.u[3],
    integrator.u[4],
    integrator.p.r0,
    integrator.p.v0,
    integrator.p.n0,
)

function termination_condition(u, t, integrator)
    r, z, vr, vz = u
    _, _, ar, az = integrator.du
    out_of_grid_condition = out_of_grid(integrator)
    failed_condition = failed(integrator, r, z)
    return out_of_grid_condition || failed_condition
end

function affect!(integrator)
    integrator.p.finished[1] = true
    terminate!(integrator)
end

function save(
    u,
    t,
    integrator,
    radiation::Radiation,
    wind::Wind,
    parameters::Parameters,
    trajectory_id,
)
    if integrator.p.finished[1]
        return 0.0
    end
    data = integrator.p.data
    r, z, vr, vz = u
    _, _, ar, az = integrator.du
    vt = sqrt(vr^2 + vz^2)
    at = sqrt(ar^2 + az^2)
    vphi = integrator.p.l0 / sqrt(r^2 + z^2)
    dvdr = at / vt
    density = compute_density(r, z, vr, vz, integrator.p)
    taux = compute_tau_xray(
        wind.density_grid,
        wind.grid_iterator,
        radiation,
        parameters,
        ri = 0.0,
        phii = 0.0,
        zi = parameters.z_xray,
        rf = r,
        zf = z,
        phif = 0.0,
    )
    ξ = compute_ionization_parameter(
        radiation.scattered_lumin_grid,
        wind.grid_iterator,
        wind.density_grid,
        r = r,
        z = z,
        number_density = density,
        tau_x = taux,
        xray_luminosity = radiation.xray_luminosity,
        Rg = radiation.bh.Rg,
        scattering_flag = NoScattering(),
        absorption_opacity = parameters.xray_opacity_flag,
        zh = parameters.z_xray,
    )
    taueff = compute_tau_eff(density, dvdr)
    forcemultiplier = compute_force_multiplier(taueff, ξ)
    push!(data[:r], r)
    push!(data[:z], z)
    push!(data[:vr], vr)
    push!(data[:vphi], vphi)
    push!(data[:vz], vz)
    push!(data[:n], density)
    push!(data[:dvdr], dvdr)
    push!(data[:taux], taux)
    push!(data[:tauuv], 0.0)
    push!(data[:xi], ξ)
    push!(data[:taueff], taueff)
    push!(data[:fm], forcemultiplier)
    return 0.0
end

function residual!(radiation::Radiation, wind, parameters, out, du, u, p, t)
    r, z, vr, vz = u
    r_dot, z_dot, vr_dot, vz_dot = du
    if r <= 0 || z < 0 # we force it to fail
        radiation_acceleration = [0.0, 0.0]
        centrifugal_term = 0.0
        gravitational_acceleration = compute_gravitational_acceleration(abs(r), abs(z))
    else
        radiation_acceleration =
            compute_radiation_acceleration(radiation, wind, parameters, du, u, p)
        centrifugal_term = p.l0^2 / r^3
        gravitational_acceleration = compute_gravitational_acceleration(r, z)
    end
    ar = gravitational_acceleration[1] + radiation_acceleration[1] + centrifugal_term
    az = gravitational_acceleration[2] + radiation_acceleration[2]
    out[1] = r_dot - vr
    out[2] = z_dot - vz
    out[3] = vr_dot - ar
    out[4] = vz_dot - az
end

function create_dae_problem(
    radiation,
    wind,
    parameters,
    residual!,
    du₀,
    u₀,
    tspan,
    integ_params,
)
    func!(out, du, u, p, t) =
        residual!(radiation, wind, parameters, out, du, u, integ_params, t)
    return DAEProblem(
        func!,
        du₀,
        u₀,
        tspan,
        integ_params,
        differential_vars = [true, true, true, true],
    )
end

function compute_radiation_acceleration(
    radiation::Radiation,
    wind::Wind,
    parameters::Parameters,
    du,
    u,
    integ_parameters::IntegratorParameters,
)
    r, z, vr, vz = u
    _, _, ar, az = du
    vt = sqrt(vr^2 + vz^2)
    at = sqrt(ar^2 + az^2)
    if vt == 0
        dvdr = 1e-6
    else
        dvdr = at / vt
    end
    density = compute_density(r, z, vr, vz, integ_parameters)
    taux = compute_tau_xray(
        wind.density_grid,
        wind.grid_iterator,
        radiation,
        parameters,
        ri = 0.0,
        phii = 0.0,
        zi = parameters.z_xray,
        rf = r,
        phif = 0.0,
        zf = z,
    )
    #taux = 0.0
    ξ = compute_ionization_parameter(
        radiation,
        wind,
        parameters,
        r = r,
        z = z,
        number_density = density,
        tau_x = taux,
    )
    #ξ = 0.0
    taueff = compute_tau_eff(density, dvdr)
    #taueff = 1.0
    forcemultiplier = compute_force_multiplier(taueff, ξ)
    disc_radiation_field = compute_disc_radiation_field(
        wind.density_grid,
        wind.grid_iterator,
        parameters,
        radiation,
        r_wind = r,
        z_wind = z,
        vr_wind = vr,
        vz_wind = vz,
    )
    force_radiation = (1.0 + forcemultiplier) * disc_radiation_field
    return force_radiation
end

function compute_initial_acceleration(
    radiation::Radiation,
    wind::Wind,
    parameters::Parameters,
    r,
    z,
    vr,
    vz,
    integ_params::IntegratorParameters,
)
    u = [r, z, vr, vz]
    du = [vr, vz, 0, 0]
    gravitational_acceleration = compute_gravitational_acceleration(r, z)
    radiation_acceleration =
        compute_radiation_acceleration(radiation, wind, parameters, du, u, integ_params)
    centrifugal_term = integ_params.l0^2 / r^3
    if r == 0
        centrifugal_term = 0.0
    end
    ar = centrifugal_term + gravitational_acceleration[1] + radiation_acceleration[1]
    az = gravitational_acceleration[2] + radiation_acceleration[2]
    # second estimation
    du = [vr, vz, ar, az]
    radiation_acceleration =
        compute_radiation_acceleration(radiation, wind, parameters, du, u, integ_params)
    ar = centrifugal_term + gravitational_acceleration[1] + radiation_acceleration[1]
    az = gravitational_acceleration[2] + radiation_acceleration[2]
    return [ar, az]
end

function compute_lines_range(
    model,
    rin,
    rfi;
    delta_mdot = 0.01,
    fill_delta = 20,
    max_tau = 25,
)
    lines_range = []
    lines_widths = []
    r_range = range(rin, rfi, length = 250)
    density_range = getn0.(Ref(model), r_range)
    r = rin
    function mass_loss_kernel(model; r)
        ridx = min(searchsortedfirst(r_range, r), length(r_range))
        n0 = density_range[ridx]
        v0 = getv0(model, r)
        return n0 * M_P * v0 * C * 2π * r * model.bh.Rg^2
    end
    function tau_kernel(model; r)
        ridx = min(searchsortedfirst(r_range, r), length(r_range))
        n0 = density_range[ridx]
        return n0 * SIGMA_T * model.bh.Rg
    end
    get_tau(delta_r) =
        quadgk(r -> tau_kernel(model, r = r), r, r + delta_r, rtol = 1e-4, atol = 0)[1]
    mass_loss(delta_r) =
        quadgk(r -> mass_loss_kernel(model, r = r), r, r + delta_r, rtol = 1e-4, atol = 0)[1]
    target_mdot = delta_mdot * compute_mass_accretion_rate(model.bh)
    function find_delta_r(model, r; delta_mdot, tau_total)
        if mass_loss(1000) < target_mdot
            delta_r_mdot = fill_delta
        else
            delta_r_mdot = find_zero(
                delta_r -> mass_loss(delta_r) - target_mdot,
                (1e-5, 1000),
                atol = 0,
                rtol = 1e-6,
            )
        end
        if tau_total < 5
            max_delta_tau = 0.05
        elseif tau_total < 10
            max_delta_tau = 0.5
        elseif tau_total < 100
            max_delta_tau = 5
        else
            max_delta_tau = 50
        end
        if get_tau(1000) < max_delta_tau
            delta_r_tau = fill_delta
        else
            delta_r_tau = find_zero(
                delta_r -> get_tau(delta_r) - max_delta_tau,
                (1e-5, 1000),
                atol = 0,
                rtol = 1e-6,
            )
        end
        delta_r = min(delta_r_tau, delta_r_mdot)
        return min(delta_r, fill_delta)
    end
    tau_total = 0.0
    while r < rfi
        delta_r = find_delta_r(model, r, delta_mdot = delta_mdot, tau_total = tau_total)
        tau_total += get_tau(delta_r)
        push!(lines_range, r + delta_r / 2)
        push!(lines_widths, delta_r)
        r += delta_r
    end
    return lines_range, lines_widths
end

compute_lines_range(model) = compute_lines_range(
    model,
    model.ic.rin,
    model.ic.rfi,
    delta_mdot = 0.01,
    fill_delta = 20,
    max_tau = 25,
)
