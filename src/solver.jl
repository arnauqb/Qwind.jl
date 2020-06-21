using DiffEqBase
using DiffEqCallbacks
using Sundials
export Parameters, initialize_solver, run_solver!, compute_density
import Qwind.out_of_grid

struct Parameters
    line::Streamline
    grid::Grid
    radiation::Radiation
    windtree::WindTree
end

function initialize_solver(
    line::Streamline,
    params::Parameters;
    atol = 1e-8,
    rtol = 1e-3,
)
    termination_callback = DiscreteCallback(
        termination_condition,
        affect!,
        save_positions = (false, false),
    )
    saving_callback =
        SavingCallback(save, SavedValues(Float64, Array{Float64,1}))
    callback_set = CallbackSet(termination_callback, saving_callback)
    a₀ = compute_initial_acceleration(
        params.radiation,
        r(line),
        z(line),
        v_r(line),
        v_z(line),
        params,
    )
    du₀ = [v_r_0(line), v_z_0(line), a₀[1], a₀[2]]
    u₀ = [r_0(line), z_0(line), v_r_0(line), v_z_0(line)]
    tspan = (0.0, 1e8)
    dae_problem =
        create_dae_problem(params.radiation, residual!, du₀, u₀, tspan, params)
    integrator =
        init(dae_problem, IDA(init_all = false), callback = callback_set)
    integrator.opts.abstol = atol
    integrator.opts.reltol = rtol
    return integrator
end

run_solver!(solver::Sundials.IDAIntegrator) = solve!(solver)
out_of_grid(params::Parameters) = out_of_grid(params.grid, params.line)


function termination_condition(u, t, integrator)
    r, z, v_r, v_z = u
    _, _, a_r, a_z = integrator.du
    out_of_grid_condition = out_of_grid(integrator.p)
    failed_condition = z < 5 && v_z < 1e-9
    return out_of_grid_condition || failed_condition
end


function affect!(integrator)
    if escaped(integrator.p.line, integrator.p.radiation.bh.M)
        print(" \U1F4A8")
    elseif failed(integrator.p.line, integrator.p.grid)
        print(" \U1F4A5")
    else
        print(" \U2753")
    end
    terminate!(integrator)
end

function save(u, t, integrator)
    r, z, v_r, v_z = u
    line = integrator.p.line
    save_position!(line, u, t)
    density = compute_density(line)
    push!(line.number_density, density)
    #quadtree_fill_timestep(integrator.p.quadtree, line)
    return [0.0]
end

function residual!(radiation::Radiation, out, du, u, p, t)
    r, z, v_r, v_z = u
    r_dot, z_dot, v_r_dot, v_z_dot = du
    if r <= 0 || z <= 0 # we force it to fail
        radiation_acceleration = [0.0, 0.0]
        centrifugal_term = 0.0
        gravitational_acceleration =
            compute_gravitational_acceleration(abs(r), abs(z), p.radiation.bh.M)
    else
        radiation_acceleration =
            compute_radiation_acceleration(p.radiation, du, u, p)
        centrifugal_term = p.line.angular_momentum^2 / r^3
        gravitational_acceleration =
            compute_gravitational_acceleration(r, z, p.radiation.bh.M)
    end
    a_r =
        gravitational_acceleration[1] +
        radiation_acceleration[1] +
        centrifugal_term
    a_z = gravitational_acceleration[2] + radiation_acceleration[2]
    out[1] = r_dot - v_r
    out[2] = z_dot - v_z
    out[3] = v_r_dot - a_r
    out[4] = v_z_dot - a_z
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
    radiation::SimpleRadiation,
    du,
    u,
    p::Parameters,
)
    r, z, v_r, v_z = u
    _, _, a_r, a_z = du
    v_t = sqrt(v_r^2 + v_z^2)
    a_t = sqrt(a_r^2 + a_z^2)
    dv_dr = a_t / v_t
    #println("r : $r z $z v_r : $v_r v_z : $v_z")
    density = compute_density(r, z, v_t, p.line)
    #tau_x = compute_xray_tau(radiation, r, z, density, r_0(p.line))
    tau_x = compute_xray_tau(
        p.windtree.quadtree,
        r,
        z,
        0,
        xray_luminosity(radiation),
        p.radiation.bh.Rg,
    )
    #println("tau x : $tau_x")
    tau_uv = compute_uv_tau(radiation, r, z, density)
    ξ = compute_ionization_parameter(radiation, r, z, density, tau_x)
    tau_eff = compute_tau_eff(density, dv_dr, p.radiation.v_th)
    force_multiplier = compute_force_multiplier(tau_eff, ξ)
    disc_radiation_field = compute_disc_radiation_field(radiation, r, z)
    force_radiation =
        (1 + force_multiplier) * exp(-tau_uv) * disc_radiation_field
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

"Updates the density of the streamline giving its current position and velocity,
using mass conservation."
function compute_density(r, z, v_t, line::Streamline)
    #@assert r >= 0
    #@assert z >= 0
    d = sqrt(r^2 + z^2)
    radial = (r_0(line) / d)^2
    v_ratio = v_0(line) / v_t
    n = number_density_0(line) * radial * v_ratio
    return n
end

compute_density(line) =
    compute_density(r(line), z(line), sqrt(v_r(line)^2 + v_z(line)^2), line)

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
    gravitational_acceleration =
        compute_gravitational_acceleration(r, z, params.radiation.bh.M)
    radiation_acceleration =
        compute_radiation_acceleration(radiation, du, u, params)
    centrifugal_term = params.line.angular_momentum^2 / r^3
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
