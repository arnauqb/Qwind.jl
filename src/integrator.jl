using DiffEqBase, DiffEqCallbacks, Sundials, Printf
using Statistics: std, mean
using Roots: find_zero, Bisection
using Distributed
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
    keys = [
        :r,
        :z,
        :vr,
        :vz,
        :n,
        #    :ar,
        #    :az,
        :fm,
        :xi,
        :dvdr,
        #    :disc_radiation_field_r,
        #    :disc_radiation_field_z,
        :xi,
        :taueff,
        :taux,
        :tauuv,
    ]
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
    return compute_density(r, z, vr, vz, parameters.r0, parameters.v0, parameters.n0)
end

function initialize_integrator(
    #radiative_transfer::RadiativeTransfer,
    #bh::BlackHole,
    #grid::Grid,
    #initial_conditions::InitialConditions,
    model,
    r0,
    linewidth;
    atol = 1e-8,
    rtol = 1e-3,
    tmax = 1e8,
    trajectory_id = -1,
    save_results = true,
)
    rt = model.rt
    bh = model.bh
    grid = model.wind_grid
    ic = model.ic
    l0 = getl0(ic, r0)
    z0 = getz0(ic, r0)
    n0 = getn0(model, r0)
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
            (u, t, integrator) -> save(u, t, integrator, rt, trajectory_id),
            saved_data,
        )
        callback_set = CallbackSet(termination_callback, saving_callback)
    else
        callback_set = CallbackSet(termination_callback)
    end
    a₀ = compute_initial_acceleration(rt, bh, r0, z0, 0, v0, params)
    du₀ = [0.0, v0, a₀[1], a₀[2]]
    u₀ = [r0, z0, 0.0, v0]
    tspan = (0.0, tmax)
    dae_problem = create_dae_problem(rt, bh, residual!, du₀, u₀, tspan, params)
    integrator = init(dae_problem, IDA(init_all = false), callback = callback_set)
    integrator.opts.abstol = atol
    integrator.opts.reltol = rtol
    return integrator
end

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

function create_and_run_integrator(
    #radiative_transfer::RadiativeTransfer,
    #bh::BlackHole,
    #grid::Grid,
    #initial_conditions::InitialConditions;
    model;
    r0,
    linewidth,
    trajectory_id,
    atol,
    rtol,
)
    integrator = initialize_integrator(
        model,
        #radiative_transfer,
        #bh,
        #grid,
        #initial_conditions,
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
    integrator.u[1] < 0.0 || integrator.u[2] < integrator.p.grid.z_min || intersects
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

function save(u, t, integrator, radiative_transfer::RadiativeTransfer, trajectory_id)
    if integrator.p.finished[1]
        return 0.0
    end
    data = integrator.p.data
    r, z, vr, vz = u
    _, _, ar, az = integrator.du
    vt = sqrt(vr^2 + vz^2)
    at = sqrt(ar^2 + az^2)
    dvdr = at / vt
    density = compute_density(r, z, vr, vz, integrator.p)
    taux = compute_xray_tau(radiative_transfer, radiative_transfer.radiation.z_xray, r, z)
    ξ = compute_ionization_parameter(radiative_transfer.radiation, r, z, density, taux)
    taueff = compute_tau_eff(density, dvdr)
    forcemultiplier = compute_force_multiplier(taueff, ξ)
    push!(data[:r], r)
    push!(data[:z], z)
    push!(data[:vr], vr)
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

function save(u, t, integrator, radiative_transfer::RERadiativeTransfer, trajectory_id)
    if integrator.p.finished[1]
        return 0.0
    end
    data = integrator.p.data
    r, z, vr, vz = u
    _, _, ar, az = integrator.du
    vt = sqrt(vr^2 + vz^2)
    at = sqrt(ar^2 + az^2)
    dvdr = at / vt
    Rg = radiative_transfer.radiation.Rg
    density = compute_density(r, z, vr, vz, integrator.p)
    taux = compute_xray_tau(radiative_transfer, r, z, density, integrator.p.r0)
    tauuv = compute_uv_tau(radiative_transfer, r, z, density, integrator.p.r0)
    ξ = compute_ionization_parameter(radiative_transfer.radiation, r, z, density, taux)
    taueff = compute_tau_eff(density, dvdr)
    forcemultiplier = compute_force_multiplier(taueff, ξ)
    disc_radiation_field = compute_disc_radiation_field(radiative_transfer, r, z)
    push!(data[:r], r)
    push!(data[:z], z)
    push!(data[:vr], vr)
    push!(data[:vz], vz)
    push!(data[:ar], ar)
    push!(data[:az], az)
    push!(data[:n], density)
    push!(data[:dvdr], dvdr)
    push!(data[:taux], taux)
    push!(data[:tauuv], tauuv)
    push!(data[:xi], ξ)
    push!(data[:taueff], taueff)
    push!(data[:fm], forcemultiplier)
    push!(data[:disc_radiation_field_r], disc_radiation_field[1])
    push!(data[:disc_radiation_field_z], disc_radiation_field[2])
    return 0.0
end

function residual!(radiative_transfer::RadiativeTransfer, bh::BlackHole, out, du, u, p, t)
    r, z, vr, vz = u
    r_dot, z_dot, vr_dot, vz_dot = du
    if r <= 0 || z < 0 # we force it to fail
        radiation_acceleration = [0.0, 0.0]
        centrifugal_term = 0.0
        gravitational_acceleration =
            compute_gravitational_acceleration(bh, abs(r), abs(z), zh = "height")
    else
        radiation_acceleration =
            compute_radiation_acceleration(radiative_transfer, du, u, p)
        centrifugal_term = p.l0^2 / r^3
        gravitational_acceleration = compute_gravitational_acceleration(bh, r, z, zh = "height")
    end
    ar = gravitational_acceleration[1] + radiation_acceleration[1] + centrifugal_term
    az = gravitational_acceleration[2] + radiation_acceleration[2]
    out[1] = r_dot - vr
    out[2] = z_dot - vz
    out[3] = vr_dot - ar
    out[4] = vz_dot - az
end

function create_dae_problem(
    radiative_transfer::RadiativeTransfer,
    bh::BlackHole,
    residual!,
    du₀,
    u₀,
    tspan,
    params,
)
    func!(out, du, u, p, t) = residual!(radiative_transfer, bh, out, du, u, p, t)
    return DAEProblem(
        func!,
        du₀,
        u₀,
        tspan,
        params,
        differential_vars = [true, true, true, true],
    )
end


function compute_radiation_acceleration(
    radiative_transfer::RERadiativeTransfer,
    du,
    u,
    p::Parameters,
)
    r, z, vr, vz = u
    _, _, ar, az = du
    vt = sqrt(vr^2 + vz^2)
    at = sqrt(ar^2 + az^2)
    dvdr = at / vt
    density = compute_density(r, z, vr, vz, p)
    taux = compute_xray_tau(radiative_transfer, r, z, density, p.r0)
    tauuv = compute_uv_tau(radiative_transfer, r, z, density, p.r0)
    ξ = compute_ionization_parameter(radiative_transfer.radiation, r, z, density, taux)
    taueff = compute_tau_eff(density, dvdr)
    forcemultiplier = compute_force_multiplier(taueff, ξ)
    disc_radiation_field = compute_disc_radiation_field(radiative_transfer, r, z)
    force_radiation = (1 + forcemultiplier) * exp(-tauuv) * disc_radiation_field
    return force_radiation
end

function compute_radiation_acceleration(
    radiative_transfer::RadiativeTransfer,
    du,
    u,
    p::Parameters,
)
    r, z, vr, vz = u
    _, _, ar, az = du
    vt = sqrt(vr^2 + vz^2)
    at = sqrt(ar^2 + az^2)
    dvdr = at / vt
    density = compute_density(r, z, vr, vz, p)
    taux = compute_xray_tau(radiative_transfer, radiative_transfer.radiation.z_xray, r, z)
    ξ = compute_ionization_parameter(radiative_transfer.radiation, r, z, density, taux)
    taueff = compute_tau_eff(density, dvdr)
    forcemultiplier = compute_force_multiplier(taueff, ξ)
    disc_radiation_field = compute_disc_radiation_field(radiative_transfer, r, z, vr, vz)
    force_radiation = (1 + forcemultiplier) * disc_radiation_field
    return force_radiation
end

function compute_initial_acceleration(
    radiative_transfer::RadiativeTransfer,
    bh::BlackHole,
    r,
    z,
    vr,
    vz,
    params::Parameters,
)
    u = [r, z, vr, vz]
    du = [vr, vz, 0, 0]
    gravitational_acceleration = compute_gravitational_acceleration(bh, r, z, zh = "height")
    radiation_acceleration =
        compute_radiation_acceleration(radiative_transfer, du, u, params)
    centrifugal_term = params.l0^2 / r^3
    ar = centrifugal_term + gravitational_acceleration[1] + radiation_acceleration[1]
    az = gravitational_acceleration[2] + radiation_acceleration[2]
    # second estimation
    du = [vr, vz, ar, az]
    radiation_acceleration =
        compute_radiation_acceleration(radiative_transfer, du, u, params)
    ar = centrifugal_term + gravitational_acceleration[1] + radiation_acceleration[1]
    az = gravitational_acceleration[2] + radiation_acceleration[2]
    return [ar, az]
end

function compute_delta_r_x(grid::DensityGrid, rc, xray_luminosity, Rg, tau_x, delta_tau)
    delta_r = find_zero(
        delta_r ->
            delta_tau - (
                compute_xray_tau(grid, 0.0, 0.0, rc + delta_r, 0.0, xray_luminosity, Rg) - tau_x
            ),
        1,
        atol = 1e-7,
        rtol = 1e-3,
    )
    return delta_r
end

function compute_delta_r_uv(grid::DensityGrid, rc, Rg, tau_uv, delta_tau)
    delta_r = find_zero(
        delta_r ->
            delta_tau - (compute_uv_tau(grid, 0.0, 0.0, rc + delta_r, 0.0, Rg) - tau_uv),
        1.0,
        atol = 1e-7,
        rtol = 1e-3,
    )
    return delta_r
end

function compute_delta_r(grid::DensityGrid, rc, xray_luminosity, Rg)
    density = get_density(grid, rc, 0.0)
    tau_x_0 = log(xray_luminosity / (1e5 * density))
    tau_x = compute_xray_tau(grid, 0.0, 0.0, rc, 0.0, xray_luminosity, Rg)
    tau_uv = compute_uv_tau(grid, 0.0, 0.0, rc, 0.0, Rg)
    delta_tau = 0.1
    delta_r = 0.0
    if tau_x < 2 * tau_x_0
        try
            delta_r = compute_delta_r_x(grid, rc, xray_luminosity, Rg, tau_x, delta_tau)
        catch
            return 0.0
        end
        return min(delta_r, 0.5)
    else
        try
            delta_r = compute_delta_r_uv(grid, rc, Rg, tau_uv, delta_tau)
            return min(delta_r, 5)
        catch
            return 0.0
        end
    end
end


function compute_lines_range_old(ic::InitialConditions, rin, rfi, Rg, xray_luminosity)
    lines_range = [rin]
    lines_widths = []
    r_range = 10 .^ range(log10(rin), log10(rfi), length = 100000)
    z_range = [0.0, 1.0]
    density_grid = zeros((length(r_range), length(z_range)))
    density_grid[:, 1] .= getn0.(Ref(ic), r_range)
    interp_grid = DensityGrid(r_range, z_range, density_grid)
    tau_x = 0.0
    tau_uv = 0.0
    fx(delta_r, rc, delta_tau, tau_x) =
        delta_tau - (
            compute_xray_tau(
                interp_grid,
                0.0,
                0.0,
                rc + delta_r,
                0.0,
                xray_luminosity,
                Rg,
            ) - tau_x
        )
    fuv(delta_r, rc, delta_tau, tau_uv) =
        delta_tau - (compute_uv_tau(interp_grid, 0.0, 0.0, rc + delta_r, 0.0, Rg) - tau_uv)
    rc = rin
    while rc < rfi
        delta_r = compute_delta_r(interp_grid, rc, xray_luminosity, Rg)
        if delta_r == 0.0
            break
        end
        push!(lines_range, rc + delta_r / 2)
        push!(lines_widths, delta_r)
        rc += delta_r
    end
    # add the same amount of trajectories for the rest
    additional_range = range(rc, rfi, step = 5) #length=length(lines_range))
    additional_widths = diff(additional_range)
    pushfirst!(additional_widths, additional_range[1] - lines_range[end])
    lines_range = vcat(lines_range, additional_range)
    lines_widths = vcat(lines_widths, additional_widths)
    # discard last radius
    lines_range = lines_range[1:(end - 1)]
    return lines_range, lines_widths
end

function compute_lines_range(model, rin, rfi, Rg, xray_luminosity)
    lines_range = [rin]
    lines_widths = []
    r_range = 10 .^ range(log10(rin), log10(rfi), length = 100000)
    z_range = [0.0, 1.0]
    density_grid = zeros((length(r_range), length(z_range)))
    density_grid[:, 1] .= getn0.(Ref(model), r_range)
    interp_grid = DensityGrid(r_range, z_range, density_grid)
    tau_x = 0.0
    tau_uv = 0.0
    fx(delta_r, rc, delta_tau, tau_x) =
        delta_tau - (
            compute_xray_tau(
                interp_grid,
                0.0,
                0.0,
                rc + delta_r,
                0.0,
                xray_luminosity,
                Rg,
            ) - tau_x
        )
    fuv(delta_r, rc, delta_tau, tau_uv) =
        delta_tau - (compute_uv_tau(interp_grid, 0.0, 0.0, rc + delta_r, 0.0, Rg) - tau_uv)
    rc = rin
    while rc < rfi
        if tau_x < 50
            tau_x = compute_xray_tau(interp_grid, 0.0, 0.0, rc, 0.0, xray_luminosity, Rg)
            if tau_x < 1
                delta_tau = 0.1
            elseif tau_x < 20
                delta_tau = 0.5
            else
                delta_tau = 1
            end
            delta_r = find_zero(
                delta_r -> fx(delta_r, rc, delta_tau, tau_x),
                0.1,
                atol = 1e-7,
                rtol = 1e-3,
            )
        else
            #elseif tau_uv < 50
            tau_uv = compute_uv_tau(interp_grid, 0.0, 0.0, rc, 0.0, Rg)
            if tau_uv < 1
                delta_tau = 0.1
            elseif tau_uv < 20
                delta_tau = 0.5
            elseif tau_uv < 50
                delta_tau = 1
            else
                delta_tau = 10
                #break
            end
            try
                delta_r = find_zero(
                    delta_r -> fuv(delta_r, rc, delta_tau, tau_uv),
                    10,
                    atol = 1e-6,
                    rtol = 1e-3,
                )
            catch
                break
            end
            #else
            #    break
        end
        delta_r = min(delta_r, 1)
        push!(lines_range, rc + delta_r / 2)
        push!(lines_widths, delta_r)
        rc += delta_r
    end
    additional_range = range(rc, rfi, step = 5)
    additional_widths = diff(additional_range)
    pushfirst!(additional_widths, additional_range[1] - lines_range[end])
    lines_range = vcat(lines_range, additional_range)
    lines_widths = vcat(lines_widths, additional_widths)
    # discard last radius
    lines_range = lines_range[1:(end - 1)]
    return lines_range, lines_widths
end


compute_lines_range(model) = compute_lines_range(
    model,
    model.ic.rin,
    model.ic.rfi,
    model.bh.Rg,
    model.rad.xray_luminosity,
)
