using DiffEqBase, DiffEqCallbacks, Sundials
using DataStructures, Printf
using Statistics: std, mean
using Distributed
import Qwind.out_of_grid, Qwind.compute_density
import Base.push!
export initialize_integrator,
    initialize_integrators,
    run_integrator!,
    run_integrators!,
    run_integrators_parallel!,
    get_initial_radii_and_linewidths,
    compute_density,
    compute_disc_radiation_field,
    compute_vt,
    compute_d,
    failed,
    Parameters,
    is_stalled,
    escaped,
    get_dense_solution_from_integrator,
    get_dense_solution_from_integrators

function make_save_data(line_id = -1)
    ret = Dict{Symbol,Any}(:line_id => line_id)
    keys = [
        :r,
        :z,
        :vr,
        :vz,
        :n,
        :ar,
        :az,
        :fm,
        :xi,
        :dvdr,
        :disc_radiation_field_r,
        :disc_radiation_field_z,
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
    line_id::Int
    r0::Float64
    z0::Float64
    v0::Float64
    n0::Float64
    l0::Float64
    lwnorm::Float64
    grid::Grid
    data::Any
end

function compute_density(r, z, vr, vz, parameters::Parameters)
    return compute_density(r, z, vr, vz, parameters.r0, parameters.v0, parameters.n0)
end

function initialize_integrator(
    radiative_transfer::RadiativeTransfer,
    grid::Grid,
    initial_conditions::InitialConditions,
    r0,
    linewidth;
    atol = 1e-8,
    rtol = 1e-3,
    tmax = 1e8,
    line_id = -1,
    save_results = true,
    #save_path = nothing,
)
    l0 = getl0(initial_conditions, r0)
    z0 = getz0(initial_conditions, r0)
    n0 = getn0(initial_conditions, r0)
    v0 = getv0(initial_conditions, r0)
    lwnorm = linewidth / r0
    termination_callback = DiscreteCallback(
        termination_condition,
        integrator -> affect!(integrator),
        save_positions = (false, false),
    )
    data = make_save_data(line_id)
    stalling_cb =
        DiscreteCallback(is_stalled, stalling_affect!, save_positions = (false, false))
    params = Parameters(line_id, r0, z0, v0, n0, l0, lwnorm, grid, data)
    if save_results
        saved_data = SavedValues(Float64, Float64)
        saving_callback = SavingCallback(
            (u, t, integrator) -> save(u, t, integrator, radiative_transfer, line_id),
            saved_data,
        )
        callback_set = CallbackSet(termination_callback, saving_callback, stalling_cb)
    else
        callback_set = CallbackSet(termination_callback, stalling_cb)
    end
    a₀ = compute_initial_acceleration(radiative_transfer, r0, z0, 0, v0, params)
    du₀ = [0.0, v0, a₀[1], a₀[2]]
    u₀ = [r0, z0, 0.0, v0]
    tspan = (0.0, tmax)
    dae_problem = create_dae_problem(radiative_transfer, residual!, du₀, u₀, tspan, params)
    integrator = init(dae_problem, IDA(init_all = false), callback = callback_set)
    integrator.opts.abstol = atol
    integrator.opts.reltol = rtol
    return integrator
end

function get_initial_radii_and_linewidths(initial_conditions::InitialConditions, Rg)
    rin = getrin(initial_conditions)
    rfi = getrfi(initial_conditions)
    nlines = getnlines(initial_conditions)
    if nlines == "auto"
        lines_range, lines_widths = compute_lines_range(initial_conditions, rin, rfi, Rg)
    else
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
    end
    return lines_range, lines_widths
end

function initialize_integrators(
    radiative_transfer::RadiativeTransfer,
    grid::Grid,
    initial_conditions::InitialConditions;
    atol = 1e-8,
    rtol = 1e-3,
    #save_path = nothing,
)
    lines_range, lines_widths = get_initial_radii_and_linewidths(
        initial_conditions,
        radiative_transfer.radiation.Rg,
    )
    #if isfile(save_path)
    #    rm(save_path)
    #end
    integrators = Array{Sundials.IDAIntegrator}(undef, length(lines_range))
    for (i, (r0, linewidth)) in enumerate(zip(lines_range, lines_widths))
        integrator = initialize_integrator(
            radiative_transfer,
            grid,
            initial_conditions,
            r0,
            linewidth,
            atol = atol,
            rtol = rtol,
            line_id = i,
            #save_path = save_path,
        )
        integrators[i] = integrator
    end
    return integrators
end

function create_and_run_integrator(
    radiative_transfer::RadiativeTransfer,
    grid::Grid,
    initial_conditions::InitialConditions;
    r0,
    linewidth,
    line_id,
    atol,
    rtol,
)
    println("Running integrator $line_id, time $(get_time())")
    flush(stdout)
    flush(stderr)
    integrator = initialize_integrator(
        radiative_transfer,
        grid,
        initial_conditions,
        r0,
        linewidth,
        atol = atol,
        rtol = rtol,
        line_id = line_id,
        #save_path = save_path,
    )
    solve!(integrator)
    println("Integrator $line_id done, time $(get_time())")
    flush(stdout)
    flush(stderr)
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

function failed(integrator::Sundials.IDAIntegrator)
    sign_changes1 = countsignchanges(integrator.p.data[:vz])
    #sign_changes2 = countsignchanges(diff(integrator.p.data[:vz]), 1e-4)
    integrator.u[1] < 0.0 ||
        integrator.u[2] < integrator.p.grid.z_min ||
        (sign_changes1 >= 2)
    #(sign_changes2 >= 2)
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
    failed_condition = failed(integrator)
    return out_of_grid_condition || failed_condition
end

function is_stalled(u, t, integrator)
    return false
    #min_length = 100
    #if length(integrator.p.data[:vz]) <= min_length
    #    return false
    #end
    #vz_std =
    #    std(integrator.p.data[:vz][(end - min_length):end]) /
    #    mean(integrator.p.data[:vz][(end - min_length):end])
    #abs(vz_std) < 0.05
end

function stalling_affect!(integrator)
    @warn "STALLING!"
    flush(stdout)
    integrator.u[2] += sign(integrator.u[4]) * 5e-2 * integrator.u[2]
    integrator.u[1] += sign(integrator.u[3]) * 5e-2 * integrator.u[1]
end

function affect!(integrator)
    if escaped(integrator)
        #println("Line $(integrator.p.line_id) escaped!")
        flush(stdout)
        flush(stderr)
        #println(" \U1F4A8")
    elseif failed(integrator)
        #println("Line $(integrator.p.line_id) failed!")
        flush(stdout)
        flush(stderr)
        #println(" \U1F4A5")
    else
        #println("Line $(integrator.p.line_id) stalled!")
        flush(stdout)
        flush(stderr)
        #println(" \U2753")
    end
    terminate!(integrator)
end

function save(u, t, integrator, radiative_transfer::RadiativeTransfer, line_id)
    data = integrator.p.data
    r, z, vr, vz = u
    _, _, ar, az = integrator.du
    vt = sqrt(vr^2 + vz^2)
    at = sqrt(ar^2 + az^2)
    dvdr = at / vt
    density = compute_density(r, z, vr, vz, integrator.p)
    taux = compute_xray_tau(radiative_transfer, r, z)
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
    push!(data[:tauuv], 0.0)
    push!(data[:xi], ξ)
    push!(data[:taueff], taueff)
    push!(data[:fm], forcemultiplier)
    push!(data[:disc_radiation_field_r], disc_radiation_field[1])
    push!(data[:disc_radiation_field_z], disc_radiation_field[2])
    return 0.0
end

function save(u, t, integrator, radiative_transfer::RERadiativeTransfer, line_id)
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

function residual!(radiative_transfer::RadiativeTransfer, out, du, u, p, t)
    r, z, vr, vz = u
    r_dot, z_dot, vr_dot, vz_dot = du
    if r <= 0 || z < 0 # we force it to fail
        radiation_acceleration = [0.0, 0.0]
        centrifugal_term = 0.0
        gravitational_acceleration = compute_gravitational_acceleration(abs(r), abs(z))
    else
        radiation_acceleration =
            compute_radiation_acceleration(radiative_transfer, du, u, p)
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
    radiative_transfer::RadiativeTransfer,
    residual!,
    du₀,
    u₀,
    tspan,
    params,
)
    func!(out, du, u, p, t) = residual!(radiative_transfer, out, du, u, p, t)
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
    taux = compute_xray_tau(radiative_transfer, r, z)
    ξ = compute_ionization_parameter(radiative_transfer.radiation, r, z, density, taux)
    taueff = compute_tau_eff(density, dvdr)
    forcemultiplier = compute_force_multiplier(taueff, ξ)
    disc_radiation_field = compute_disc_radiation_field(radiative_transfer, r, z)
    fc = flux_correction(radiative_transfer.radiation.flux_correction, vt)
    force_radiation = (1 + forcemultiplier) * disc_radiation_field * fc
    return force_radiation
end

function compute_initial_acceleration(
    radiative_transfer::RadiativeTransfer,
    r,
    z,
    vr,
    vz,
    params::Parameters,
)
    u = [r, z, vr, vz]
    du = [vr, vz, 0, 0]
    gravitational_acceleration = compute_gravitational_acceleration(r, z)
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

function get_dense_solution_from_integrator(integrator, n_timesteps = 10000)
    integration_time = unique(integrator.sol.t)
    tmin = integration_time[2]
    tmax = integration_time[end - 1]
    if tmax <= tmin
        return [], [], [], [], [], []
    end
    t_range = 10 .^ range(log10(tmin), log10(tmax), length = n_timesteps)
    i = 1
    while (t_range[end] > integrator.sol.t[end]) && length(t_range) > 0
        t_range = t_range[1:(end - i)]
        i += 1
    end
    if length(t_range) == 0
        return [], [], [], [], [], []
    end
    dense_solution = integrator.sol(t_range)
    r_dense = dense_solution[1, :]
    z_dense = dense_solution[2, :]
    line_width_dense = r_dense .* integrator.p.lwnorm
    density_dense =
        compute_density.(
            dense_solution[1, :],
            dense_solution[2, :],
            dense_solution[3, :],
            dense_solution[4, :],
            Ref(integrator.p),
        )
    zmax_dense = maximum(dense_solution[2, :]) .* ones(length(dense_solution[2, :]))
    z0_dense = dense_solution[2, 1] .* ones(length(dense_solution[2, :]))
    return r_dense, z_dense, zmax_dense, z0_dense, line_width_dense, density_dense
end

function get_dense_solution_from_integrators(integrators, n_timesteps = 10000)
    r_dense = Float64[]
    z_dense = Float64[]
    zmax_dense = Float64[]
    z0_dense = Float64[]
    line_width_dense = Float64[]
    density_dense = Float64[]
    @info "Getting dense solution from integrators..."
    flush(stdout)
    for integrator in integrators
        rp, zp, zmaxp, z0p, lwp, densityp =
            get_dense_solution_from_integrator(integrator, n_timesteps)
        r_dense = vcat(r_dense, rp)
        z_dense = vcat(z_dense, zp)
        zmax_dense = vcat(zmax_dense, zmaxp)
        z0_dense = vcat(z0_dense, z0p)
        line_width_dense = vcat(line_width_dense, lwp)
        density_dense = vcat(density_dense, densityp)
    end
    @info "Done"
    flush(stdout)
    return r_dense, z_dense, zmax_dense, z0_dense, line_width_dense, density_dense
end

function compute_lines_range(ic::InitialConditions, rin, rfi, Rg)
    rc = rin
    lines_range = []
    lines_widths = []
    tau_x = 0
    tau_uv = 0
    while rc < rfi
        if tau_x < 1
            delta_r = 0.1 / (SIGMA_T * Rg * 100 * getn0(ic, rc))
            tau_x += 0.1
            tau_uv += 0.1 / 100
        elseif (tau_x > 1) && (tau_uv < 1)
            delta_r = 0.1 / (SIGMA_T * Rg * getn0(ic, rc))
            tau_x += 10
            tau_uv += 0.1
        else
            delta_r = 1 / (SIGMA_T * Rg * getn0(ic, rc))
        end
        push!(lines_range, rc + delta_r / 2)
        push!(lines_widths, delta_r)
        rc += delta_r
    end
    return lines_range, lines_widths
end
