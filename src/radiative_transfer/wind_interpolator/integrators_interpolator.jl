import Interpolations: interpolate
export reduce_integrators

"""
Reduces line to ensure it's monothonically increasing in z.
This allows to interpolate as a function of z.
"""
function reduce_line(r::T, z::T, vr::T, vz::T, n::T) where {T<:Vector{<:AbstractFloat}}
    rs = [r[1]]
    zs = [z[1]]
    vrs = [vr[1]]
    vzs = [vz[1]]
    ns = [n[1]]
    for (rp, zp, vrp, vzp, np) in zip(r, z, vr, vz, n)
        if (zp > zs[end])
            push!(rs, rp)
            push!(zs, zp)
            push!(vrs, vrp)
            push!(vzs, vzp)
            push!(ns, np)
        end
    end
    return rs, zs, vrs, vzs, ns
end

function reduce_line(integrator::Sundials.IDAIntegrator; n_timesteps = 10000, log = true)
    dense_integrator = DenseIntegrator(integrator, n_timesteps = n_timesteps, log = log)
    return reduce_line(
        dense_integrator.r,
        dense_integrator.z,
        dense_integrator.vr,
        dense_integrator.vz,
        dense_integrator.n,
    )
end

struct DenseLine{T<:AbstractFloat}
    r::Vector{T}
    z::Vector{T}
    vr::Vector{T}
    vz::Vector{T}
    n::Vector{T}
    zmin::T
    zmax::T
    interpolator::Any
    function DenseLine(integrator::Sundials.IDAIntegrator; n_timesteps = 10000, log=true)
        r, z, vr, vz, n = reduce_line(integrator, n_timesteps = n_timesteps, log=log)
        interpolator = Interpolations.interpolate((z,), r, Gridded(Linear()))
        return new{typeof(r[1])}(r, z, vr, vz, n, minimum(z), maximum(z), interpolator)
    end
end

function interpolate(line::DenseLine, z)
    if z < line.zmin
        return line.zmin
    elseif z > line.zmax
        return line.zmax
    else
        return line.interpolator(z)
    end
end

function find_intersection(integ1::DenseLine, integ2::DenseLine)
    zmin = max(integ1.zmin, integ2.zmin)
    zmax = min(integ1.zmax, integ2.zmax)
    (zmax < zmin) && return false
    z_range = range(zmin, zmax, length = 2000)
    f(z) = interpolate(integ1, z) - interpolate(integ2, z)
    previous_sign = sign(f(z_range[1]))
    i = 1
    for (i, z) in enumerate(z_range)
        current_sign = sign(f(z))
        if current_sign != previous_sign
            return z_range[i - 1]
        end
        previous_sign = current_sign
    end
    return Inf
end

function reduce_integrators(integrators; n_timesteps = 10000, log = true)
    rs = Float64[]
    zs = Float64[]
    vrs = Float64[]
    vzs = Float64[]
    ns = Float64[]
    @info "Reducing integrators..."
    flush()
    dense_lines = DenseLine.(integrators, n_timesteps=n_timesteps, log=log)
    for (i, line) in enumerate(dense_lines)
        zint = Inf
        for (j, line2) in enumerate(dense_lines)
            if j <= i
                continue
            end
            zint = min(zint, find_intersection(line, line2))
        end
        z_idx = searchsorted_nearest(line.z, zint)
        rs = vcat(rs, line.r[1:z_idx])
        zs = vcat(zs, line.z[1:z_idx])
        vrs = vcat(vrs, line.vr[1:z_idx])
        vzs = vcat(vzs, line.vz[1:z_idx])
        ns = vcat(ns, line.n[1:z_idx])
    end
    @info "Done"
    flush()
    return rs, zs, vrs, vzs, ns
end

function reduce_integrators_old(integrators; n_timesteps = 10000, log = true)
    rs = Float64[]
    zs = Float64[]
    vrs = Float64[]
    vzs = Float64[]
    ns = Float64[]
    for integ in integrators
        r, z, vr, vz, n = reduce_line(integ, n_timesteps = n_timesteps, log = true)
        rs = vcat(rs, r)
        zs = vcat(zs, z)
        vrs = vcat(vrs, vr)
        vzs = vcat(vzs, vz)
        ns = vcat(ns, n)
    end
    return rs, zs, vrs, vzs, ns
end


