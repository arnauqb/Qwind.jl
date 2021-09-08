"""
To obtain a streamline, we usually take the integration result, and resolve the intersections
with other trajectories. A streamline thus is a trajectory that does not intersect any other (including itself).
"""

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
Base.getindex(sl::Streamlines, i) = s1.streamlines[i]
Base.length(sl::Streamlines) = length(sl.streamlines)

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
        sqrt.(integrator.p.data[:r].^2 + integrator.p.data[:z] .^2) .* integ.p.lwnorm
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
