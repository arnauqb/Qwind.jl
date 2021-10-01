"""
To obtain a streamline, we usually take the integration result, and resolve the intersections
with other trajectories. A streamline thus is a trajectory that does not intersect any other (including itself).
"""

using Sundials
export Streamline, Streamlines

struct Streamline{T<:Vector{<:AbstractFloat}}
    id::Int
    t::T
    r::T
    z::T
    vr::T
    vphi::T
    vz::T
    n::T
    width::T
end

struct Streamlines
    streamlines::Vector{Streamline}
end
Base.getindex(sl::Streamlines, i) = sl.streamlines[i]
Base.length(sl::Streamlines) = length(sl.streamlines)
Base.iterate(sl::Streamlines) = iterate(sl.streamlines)
Base.iterate(sl::Streamlines, i) = iterate(sl.streamlines, i)

function Streamline(integrator::Sundials.IDAIntegrator)
    return Streamline(
        integrator.p.id,
        integrator.sol.t[1:(end - 2)],
        integrator.p.data[:r],
        integrator.p.data[:z],
        integrator.p.data[:vr],
        integrator.p.data[:vphi],
        integrator.p.data[:vz],
        integrator.p.data[:n],
        sqrt.(integrator.p.data[:r].^2 + integrator.p.data[:z] .^2) .* integrator.p.lwnorm
    )
end

function Streamline(tdata::Dict)
    return Streamline(
        tdata["id"],
        tdata["t"],
        tdata["r"],
        tdata["z"],
        tdata["vr"],
        tdata["vphi"],
        tdata["vz"],
        tdata["n"],
        tdata["line_width"]
    )
end

function Streamlines(tsdata::Dict)
    t_ids = sort(parse.(Int, keys(tsdata)))
    t_ids = string.(t_ids)
    ret = []
    for i in t_ids
        tdata = tsdata[i]
        streamline = Streamline(tdata)
        push!(ret, streamline)
    end
    return Streamlines(ret) 
end

function Streamlines(h5_path::String, it_num)
    it_name = @sprintf "iteration_%03d" it_num
    tsdata = h5open(h5_path, "r") do file
        read(file, it_name * "/streamlines")
    end
    return Streamlines(tsdata)
end

function Streamlines(h5_path::String)
    it_keys = h5open(h5_path, "r") do file
        keys(read(file))
    end
    it_nums = [parse(Int, split(key, "_")[end]) for key in it_keys]
    return Streamlines(h5_path, maximum(it_nums))
end

function slice_streamline(streamline::Streamline; in = 1, fi = nothing)
    (fi === nothing) && (fi = length(streamline.t))
    return Streamline(
        streamline.id,
        streamline.t[in:fi],
        streamline.r[in:fi],
        streamline.z[in:fi],
        streamline.vr[in:fi],
        streamline.vphi[in:fi],
        streamline.vz[in:fi],
        streamline.n[in:fi],
    )
end

function escaped(streamline::Streamline)
    vt = sqrt.(streamline.vr .^2 + streamline.vz .^2 )
    d = sqrt.(streamline.r .^2 + streamline.z .^ 2)
    return maximum(vt .> compute_escape_velocity.(d))
end

function reduce_streamlines(streamlines::Streamlines)
    r = reduce(vcat, [streamline.r for streamline in streamlines])
    z = reduce(vcat, [streamline.z for streamline in streamlines])
    vr = reduce(vcat, [streamline.vr for streamline in streamlines])
    vphi = reduce(vcat, [streamline.vphi for streamline in streamlines])
    vz = reduce(vcat, [streamline.vz for streamline in streamlines])
    n = reduce(vcat, [streamline.n for streamline in streamlines])
    return r, z, vr, vphi, vz, n
end

function interpolate_integrator(
    integrator::Sundials.IDAIntegrator;
    n_timesteps::Int = 1000,
    max_time = nothing,
    log = true,
)
    integration_time = unique(integrator.sol.t)
    tmin = integration_time[2]
    if max_time == nothing
        max_time = integration_time[end - 1]
    end
    if max_time <= tmin
        return new{Vector{Float64}}([0.0], [0.0], [0.0], [0.0], [0.0])
    end
    if log
        t_range = 10 .^ range(log10(tmin), log10(max_time), length = n_timesteps)
    else
        t_range = collect(range(tmin, max_time, length = n_timesteps))
    end
    i = 1
    while (t_range[end] > integrator.sol.t[end]) && length(t_range) > 0
        t_range = t_range[1:(end - i)]
        i += 1
    end
    if length(t_range) == 0
        return new{Vector{Float64}}([0.0], [0.0], [0.0], [0.0], [0.0])
    end
    dense_solution = integrator.sol(t_range)
    r_dense = dense_solution[1, :]
    z_dense = dense_solution[2, :]
    vr_dense = dense_solution[3, :]
    vz_dense = dense_solution[4, :]
    line_width_dense = sqrt.(r_dense .^ 2 + z_dense .^2) .* integrator.p.lwnorm
    vphi_dense = integrator.p.l0 ./ r_dense
    density_dense =
        compute_density.(r_dense, z_dense, vr_dense, vz_dense, Ref(integrator.p))
    return Streamline(
        integrator.p.id,
        t_range,
        r_dense,
        z_dense,
        vr_dense,
        vphi_dense,
        vz_dense,
        density_dense,
        line_width_dense,
    )
end

function interpolate_integrators(
    integrators::Vector{<:Sundials.IDAIntegrator};
    n_timesteps::Int = 1000,
    max_times = nothing,
    log = true,
)
    ret = []
    for integrator in integrators
        max_time = max_times[integrator.p.id]
        push!(
            ret,
            interpolate_integrator(
                integrator,
                n_timesteps = n_timesteps,
                max_time = max_time,
                log = log,
            ),
        )
    end
    return Streamlines(ret)
end
