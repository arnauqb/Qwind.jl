using DiffEqBase, DiffEqCallbacks, Sundials, Printf
using Statistics: std, mean
using Roots: find_zero, Bisection
using Distributed
using QuadGK
import Qwind.out_of_grid, Qwind.compute_density
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
    Parameters,
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

struct Parameters
    id::Int
    r0::Float64
    z0::Float64
    v0::Float64
    n0::Float64
    l0::Float64
    lwnorm::Float64
    grid::Grid
    data::Any
    finished::Vector{Bool}
end

function compute_density(r, z, vr, vz, parameters::Parameters)
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
    grid,
    ic,
    r0,
    linewidth;
    atol = 1e-8,
    rtol = 1e-3,
    tmax = 1e8,
    trajectory_id = -1,
    save_results = true,
)
    l0 = getl0(ic, r0)
    z0 = getz0(ic, r0)
    n0 = getn0(ic, rad, r0)
    v0 = getv0(ic, r0)
    lwnorm = linewidth / r0
    termination_callback = DiscreteCallback(
        termination_condition,
        integrator -> affect!(integrator),
        save_positions = (false, false),
    )
    data = make_save_data(trajectory_id)
    params = Parameters(trajectory_id, r0, z0, v0, n0, l0, lwnorm, grid, data, [false])
    if save_results
        saved_data = SavedValues(Float64, Float64)
        saving_callback = SavingCallback(
            (u, t, integrator) -> save(u, t, integrator, rad, trajectory_id),
            saved_data,
        )
        callback_set = CallbackSet(termination_callback, saving_callback)
    else
        callback_set = CallbackSet(termination_callback)
    end
    a₀ = compute_initial_acceleration(rad, r0, z0, 0, v0, params)
    du₀ = [0.0, v0, a₀[1], a₀[2]]
    u₀ = [r0, z0, 0.0, v0]
    tspan = (0.0, tmax)
    dae_problem = create_dae_problem(rad, residual!, du₀, u₀, tspan, params)
    integrator = init(dae_problem, IDA(init_all = false), callback = callback_set)
    integrator.opts.abstol = atol
    integrator.opts.reltol = rtol
    return integrator
end

initialize_integrator(
    model,
    r0,
    linewidth;
    atol = 1e-8,
    rtol = 1e-3,
    tmax = 1e8,
    trajectory_id = -1,
    save_results = true,
) = initialize_integrator(
    model.rad,
    model.wind_grid,
    model.ic,
    r0,
    linewidth,
    atol = atol,
    rtol = rtol,
    tmax = tmax,
    trajectory_id = trajectory_id,
    save_results = save_results,
)

function get_initial_radii_and_linewidths(
    initial_conditions::InitialConditions,
    xray_luminosity,
    Rg,
)
    rin = getrin(initial_conditions)
    rfi = getrfi(initial_conditions)
    nlines = getnlines(initial_conditions)
    if initial_conditions.logspaced
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

function create_and_run_integrator(model; r0, linewidth, trajectory_id, atol, rtol)
    integrator = initialize_integrator(
        model,
        r0,
        linewidth,
        atol = atol,
        rtol = rtol,
        trajectory_id = trajectory_id,
    )
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
    out_of_grid(integrator.p.grid, integrator.u[1], integrator.u[2])

function escaped(integrator::Sundials.IDAIntegrator)
    compute_vt(integrator) > compute_escape_velocity(compute_d(integrator))
end

function failed(integrator::Sundials.IDAIntegrator, r, z)
    intersects = self_intersects(integrator, r, z)
    too_long = length(integrator.p.data[:r]) > 500
    integrator.u[1] < 0.0 ||
        integrator.u[2] < max(integrator.p.grid.z_min, integrator.p.z0) ||
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

function save(u, t, integrator, radiation::Radiation, trajectory_id)
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
    taux = compute_tau_xray(radiation, r = r, z = z)
    ξ = compute_ionization_parameter(radiation, r, z, vr, vz, density, taux)
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

function residual!(radiation::Radiation, out, du, u, p, t)
    r, z, vr, vz = u
    r_dot, z_dot, vr_dot, vz_dot = du
    if r <= 0 || z < 0 # we force it to fail
        radiation_acceleration = [0.0, 0.0]
        centrifugal_term = 0.0
        gravitational_acceleration = compute_gravitational_acceleration(abs(r), abs(z))
    else
        radiation_acceleration = compute_radiation_acceleration(radiation, du, u, p)
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

function create_dae_problem(radiation, residual!, du₀, u₀, tspan, params)
    func!(out, du, u, p, t) = residual!(radiation, out, du, u, p, t)
    return DAEProblem(
        func!,
        du₀,
        u₀,
        tspan,
        params,
        differential_vars = [true, true, true, true],
    )
end

function compute_radiation_acceleration(radiation::Radiation, du, u, p::Parameters)
    r, z, vr, vz = u
    _, _, ar, az = du
    vt = sqrt(vr^2 + vz^2)
    at = sqrt(ar^2 + az^2)
    dvdr = at / vt
    density = compute_density(r, z, vr, vz, p)
    taux = compute_tau_xray(radiation, r = r, z = z)
    ξ = compute_ionization_parameter(radiation, r, z, vr, vz, density, taux)
    taueff = compute_tau_eff(density, dvdr)
    forcemultiplier = compute_force_multiplier(taueff, ξ)
    disc_radiation_field = compute_disc_radiation_field(
        radiation,
        r = r,
        z = z,
        vr = vr,
        vz = vz,
        rtol = radiation.disk_integral_rtol,
    )
    force_radiation = (1 + forcemultiplier) * disc_radiation_field
    return force_radiation
end

function compute_initial_acceleration(
    radiation::Radiation,
    r,
    z,
    vr,
    vz,
    params::Parameters,
)
    u = [r, z, vr, vz]
    du = [vr, vz, 0, 0]
    gravitational_acceleration = compute_gravitational_acceleration(r, z)
    radiation_acceleration = compute_radiation_acceleration(radiation, du, u, params)
    centrifugal_term = params.l0^2 / r^3
    if r == 0
        centrifugal_term = 0.0
    end
    ar = centrifugal_term + gravitational_acceleration[1] + radiation_acceleration[1]
    az = gravitational_acceleration[2] + radiation_acceleration[2]
    # second estimation
    du = [vr, vz, ar, az]
    radiation_acceleration = compute_radiation_acceleration(radiation, du, u, params)
    ar = centrifugal_term + gravitational_acceleration[1] + radiation_acceleration[1]
    az = gravitational_acceleration[2] + radiation_acceleration[2]
    return [ar, az]
end

function compute_lines_range(model, rin, rfi; delta_mdot=0.01, fill_delta = 20, max_tau=0.1)
    lines_range = []
    lines_widths = []
    r = rin
    function mass_loss_kernel(model; r)
        n0 = getn0(model, r)
        v0 = getv0(model, r)
        return n0 * M_P * v0 * C * 2π * r * model.bh.Rg^2
    end
    function tau_kernel(model; r)
        n0 = getn0(model, r)
        return n0 * SIGMA_T * model.bh.Rg
    end
    get_tau(delta_r) = quadgk(r -> tau_kernel(model, r=r), r, r+delta_r, rtol=1e-2, atol=0)[1]
    mass_loss(delta_r) = quadgk(r-> mass_loss_kernel(model, r=r), r, r + delta_r, rtol=1e-2, atol=0)[1]
    function find_delta_r(model, r; delta_mdot, tau_total)
        if tau_total > 10
            target = delta_mdot * compute_mass_accretion_rate(model.bh)
            if mass_loss(1000) < target
                return fill_delta
            end
            delta_r = find_zero(delta_r -> mass_loss(delta_r) - target, (1e-5, 1000), atol=0, rtol=1e-4)
        else
            target = max_tau
            delta_r = find_zero(delta_r -> get_tau(delta_r) - target, (1e-5, 1000), atol=0, rtol=1e-4)
        end
        return delta_r
    end
    tau_total = 0.0
    while r < rfi
        delta_r = find_delta_r(model, r, delta_mdot=delta_mdot, tau_total=tau_total)
        tau_total += get_tau(delta_r)
        push!(lines_range, r + delta_r / 2)
        push!(lines_widths, delta_r)
        r += delta_r
    end
    return lines_range, lines_widths
end

#function compute_lines_range(model, rin, rfi, Rg, xray_luminosity)
#    lines_range = [rin]
#    lines_widths = []
#    r_range = 10 .^ range(log10(rin), log10(rfi), length = 100000)
#    z_range = [0.0, 1.0]
#    density_grid = zeros((length(r_range), length(z_range)))
#    density_grid[:, 1] .= getn0.(Ref(model), r_range)
#    interp_grid = DensityGrid(r_range, z_range, density_grid)
#    tau_x = 0.0
#    tau_uv = 0.0
#    xr_opacity = model.rad.xray_opacity
#    fx(delta_r, rc, delta_tau, tau_x) =
#        delta_tau - (
#            compute_tau_xray(
#                interp_grid,
#                xr_opacity,
#                ri = 1e-6,
#                zi = 0.0,
#                rf = rc + delta_r,
#                zf = 0.0,
#                xray_luminosity = xray_luminosity,
#                Rg = Rg,
#            ) - tau_x
#        )
#    fuv(delta_r, rc, delta_tau, tau_uv) =
#        delta_tau - (
#            compute_tau_uv(
#                interp_grid,
#                ri = 0.0,
#                phii = 0.0,
#                zi = 0.0,
#                rf = rc + delta_r,
#                phif = 0.0,
#                zf = 0.0,
#                Rg = Rg,
#            ) - tau_uv
#        )
#    rc = rin
#    while rc < rfi
#        if tau_x < 50
#            tau_x = compute_tau_xray(
#                interp_grid,
#                xr_opacity,
#                ri = 0.0,
#                zi = 0.0,
#                rf = rc,
#                zf = 0.0,
#                xray_luminosity = xray_luminosity,
#                Rg = Rg,
#            )
#            if tau_x < 5
#                delta_tau = 0.1
#            elseif tau_x < 20
#                delta_tau = 0.5
#            else
#                delta_tau = 1
#            end
#            try
#                delta_r = find_zero(
#                    delta_r -> fx(delta_r, rc, delta_tau, tau_x),
#                    0.1,
#                    atol = 1e-7,
#                    rtol = 1e-3,
#                )
#            catch
#                break
#            end
#        else
#            tau_uv = compute_tau_uv(
#                interp_grid,
#                ri = 0.0,
#                phii = 0.0,
#                zi = 0.0,
#                rf = rc,
#                phif = 0.0,
#                zf = 0.0,
#                Rg = Rg,
#            )
#            if tau_uv < 5
#                delta_tau = 0.1
#            elseif tau_uv < 20
#                delta_tau = 0.5
#            elseif tau_uv < 50
#                delta_tau = 1
#            elseif tau_uv < 1000
#                delta_tau = 10
#            else
#                break
#            end
#            try
#                delta_r = find_zero(
#                    delta_r -> fuv(delta_r, rc, delta_tau, tau_uv),
#                    (0, 2),
#                    atol = 1e-6,
#                    rtol = 1e-3,
#                )
#            catch
#                break
#            end
#            if delta_r < 0
#                break
#            end
#            #else
#            #    break
#        end
#        delta_r = min(delta_r, 1)
#        push!(lines_range, rc + delta_r / 2)
#        push!(lines_widths, delta_r)
#        rc += delta_r
#    end
#    if rc < 200
#        additional_range = range(rc, 200.0, step = 2)
#        additional_range = vcat(additional_range, range(201.0, rfi, step = 5))
#    else
#        additional_range = range(rc, rfi, step = 5)
#    end
#    additional_widths = diff(additional_range)
#    pushfirst!(additional_widths, additional_range[1] - lines_range[end])
#    lines_range = vcat(lines_range, additional_range)
#    lines_widths = vcat(lines_widths, additional_widths)
#    # discard last radius
#    lines_range = lines_range[1:(end - 1)]
#    return lines_range, lines_widths
#end


compute_lines_range(model) = compute_lines_range(
    model,
    model.ic.rin,
    model.ic.rfi,
    delta_mdot=0.01,
    max_tau=0.1,
    fill_delta=20,
)
