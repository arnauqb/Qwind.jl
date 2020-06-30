using DiffEqBase, DiffEqCallbacks, Sundials
using DataStructures
import Qwind.out_of_grid, Qwind.compute_density
export initialize_integrator,
    initialize_integrators,
    run_integrator!,
    run_integrators!,
    compute_density,
    compute_disc_radiation_field,
    compute_vt,
    compute_d,
    failed,
    Parameters,
    SavedData

struct SavedData
    r::Vector{Float64}
    z::Vector{Float64}
    vr::Vector{Float64}
    vz::Vector{Float64}
    n::Vector{Float64}
    ar::Vector{Float64}
    az::Vector{Float64}
    fm::Vector{Float64}
    xi::Vector{Float64}
end

SavedData() = SavedData(
    Float64[],
    Float64[],
    Float64[],
    Float64[],
    Float64[],
    Float64[],
    Float64[],
    Float64[],
    Float64[],
)

struct Parameters
    r0::Float64
    z0::Float64
    v0::Float64
    n0::Float64
    l0::Float64
    lwnorm::Float64
    grid::Grid
    data::SavedData
end


function compute_density(r, z, vr, vz, parameters::Parameters)
    return compute_density(
        r,
        z,
        vr,
        vz,
        parameters.r0,
        parameters.v0,
        parameters.n0,
    )
end

function initialize_integrator(
    radiation::Radiation,
    grid::Grid,
    initial_conditions::InitialConditions,
    r0,
    linewidth;
    atol = 1e-8,
    rtol = 1e-3,
    tmax = 1e8
)
    l0 = getl0(initial_conditions, r0)
    z0 = getz0(initial_conditions, r0)
    n0 = getn0(initial_conditions, r0)
    v0 = getv0(initial_conditions, r0)
    lwnorm = linewidth / r0
    termination_callback = DiscreteCallback(
        termination_condition,
        affect!,
        save_positions = (false, false),
    )
    data = SavedData()
    saving_callback = SavingCallback(save, SavedValues(Float64, Float64))
    params = Parameters(r0, z0, v0, n0, l0, lwnorm, grid, data)
    callback_set = CallbackSet(termination_callback, saving_callback)
    a₀ = compute_initial_acceleration(radiation, r0, z0, 0, v0, params)
    du₀ = [0.0, v0, a₀[1], a₀[2]]
    u₀ = [r0, z0, 0.0, v0]
    tspan = (0.0, tmax)
    dae_problem =
        create_dae_problem(radiation, residual!, du₀, u₀, tspan, params)
    integrator =
        init(dae_problem, IDA(init_all = false), callback = callback_set)
    integrator.opts.abstol = atol
    integrator.opts.reltol = rtol
    return integrator
end

function initialize_integrators(
    radiation::Radiation,
    grid::Grid,
    initial_conditions::InitialConditions;
    atol = 1e-8,
    rtol = 1e-3,
)
    rin = getrin(initial_conditions)
    rfi = getrfi(initial_conditions)
    nlines = getnlines(initial_conditions)
    if initial_conditions.logspaced
        linedelimiters =
            10 .^ range(log10(rin), log10(rfi), length = nlines + 1)
        linesrange = []
        for i = 1:nlines
            r0 =
                linedelimiters[i] +
                (linedelimiters[i+1] - linedelimiters[i]) / 2.0
            push!(linesrange, r0)
        end
        lineswidths = diff(linedelimiters)
    else
        dr = (rfi - rin) / nlines
        linesrange = [rin + (i + 0.5) * dr for i = 0:nlines-1]
        lineswidths = diff([linesrange; rfi + dr / 2])
    end
    integrators = []
    for (i, (r0, linewidth)) in enumerate(zip(linesrange, lineswidths))
        integrator = initialize_integrator(
            radiation,
            grid,
            initial_conditions,
            r0,
            linewidth,
            atol = atol,
            rtol = rtol,
        )
        push!(integrators, integrator)
    end
    return integrators
end

run_integrator!(integrator::Sundials.IDAIntegrator) = solve!(integrator)
function run_integrators!(integrators::Vector)
    for integrator in integrators
        run_integrator!(integrator)
    end
end

compute_vt(integrator::Sundials.IDAIntegrator) =
    sqrt(integrator.u[3]^2 + integrator.u[4]^2)

compute_d(integrator::Sundials.IDAIntegrator) =
    sqrt(integrator.u[1]^2 + integrator.u[2]^2)

out_of_grid(integrator::Sundials.IDAIntegrator) =
    out_of_grid(integrator.p.grid, integrator.u[1], integrator.u[2])

function escaped(integrator::Sundials.IDAIntegrator)
    compute_vt(integrator) > compute_escape_velocity(compute_d(integrator))
end

function failed(integrator::Sundials.IDAIntegrator)
    integrator.u[1] < 0.0 || integrator.u[2] < integrator.p.grid.z_min
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

function affect!(integrator)
    if escaped(integrator)
        print(" \U1F4A8")
    elseif failed(integrator)
        print(" \U1F4A5")
    else
        print(" \U2753")
    end
    terminate!(integrator)
end

function save(u, t, integrator)
    data = integrator.p.data
    push!(data.r, u[1])
    push!(data.z, u[2])
    push!(data.vr, u[3])
    push!(data.vz, u[4])
    push!(data.n, compute_density(integrator))
    return 0.0
end

function residual!(radiation::Radiation, out, du, u, p, t)
    r, z, vr, vz = u
    r_dot, z_dot, vr_dot, vz_dot = du
    if r <= 0 || z < 0 # we force it to fail
        radiation_acceleration = [0.0, 0.0]
        centrifugal_term = 0.0
        gravitational_acceleration =
            compute_gravitational_acceleration(abs(r), abs(z))
    else
        radiation_acceleration =
            compute_radiation_acceleration(radiation, du, u, p)
        centrifugal_term = p.l0^2 / r^3
        gravitational_acceleration = compute_gravitational_acceleration(r, z)
    end
    ar =
        gravitational_acceleration[1] +
        radiation_acceleration[1] +
        centrifugal_term
    az = gravitational_acceleration[2] + radiation_acceleration[2]
    out[1] = r_dot - vr
    out[2] = z_dot - vz
    out[3] = vr_dot - ar
    out[4] = vz_dot - az
end

function create_dae_problem(
    radiation::Radiation,
    residual!,
    du₀,
    u₀,
    tspan,
    params,
)
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


function compute_radiation_acceleration(
    radiation::RERadiation,
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
    taux = compute_xray_tau(radiation, r, z, density, p.r0)
    tauuv = compute_uv_tau(radiation, r, z, density, p.r0)
    ξ = compute_ionization_parameter(radiation, r, z, density, taux)
    taueff = compute_tau_eff(density, dvdr)
    forcemultiplier = compute_force_multiplier(taueff, ξ)
    disc_radiation_field = compute_disc_radiation_field(
        radiation,
        r,
        z,
        radiation.bh.efficiency,
        radiation.bh.isco,
    )
    force_radiation = (1 + forcemultiplier) * exp(-tauuv) * disc_radiation_field
    return force_radiation
end

function compute_radiation_acceleration(
    radiation::QsosedRadiation,
    du,
    u,
    p::Parameters,
)
    r, z, v_r, v_z = u
    _, _, a_r, a_z = du
    v_t = sqrt(v_r^2 + v_z^2)
    a_t = sqrt(a_r^2 + a_z^2)
    dv_dr = a_t / v_t
    density = compute_density(r, z, v_t, p.line)
    tau_x = compute_xray_tau(radiation, r, z, density)
    ξ = compute_ionization_parameter(radiation, r, z, density, tau_x)
    tau_eff = compute_tau_eff(density, dv_dr, p.radiation.v_th)
    force_multiplier = compute_force_multiplier(tau_eff, ξ)
    disc_radiation_field =
        compute_disc_radiation_field(radiation, r, z, radiation.bh.M)
    force_radiation = (1 + force_multiplier) * disc_radiation_field
    return force_radiation
end

function compute_initial_acceleration(
    radiation::Radiation,
    r,
    z,
    v_r,
    v_z,
    params::Parameters,
)
    u = [r, z, v_r, v_z]
    du = [v_r, v_z, 0, 0]
    gravitational_acceleration = compute_gravitational_acceleration(r, z)
    radiation_acceleration =
        compute_radiation_acceleration(radiation, du, u, params)
    centrifugal_term = params.l0^2 / r^3
    a_r =
        centrifugal_term +
        gravitational_acceleration[1] +
        radiation_acceleration[1]
    a_z = gravitational_acceleration[2] + radiation_acceleration[2]
    # second estimation
    du = [v_r, v_z, a_r, a_z]
    radiation_acceleration =
        compute_radiation_acceleration(radiation, du, u, params)
    a_r =
        centrifugal_term +
        gravitational_acceleration[1] +
        radiation_acceleration[1]
    a_z = gravitational_acceleration[2] + radiation_acceleration[2]
    return [a_r, a_z]
end
