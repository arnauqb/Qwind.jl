import Base: intersect!, intersect
using LinearAlgebra, ProgressMeter
export Segment, intersect, reduce_integrators, self_intersects

struct Segment{T<:AbstractFloat}
    p1x::T
    p1y::T
    p2x::T
    p2y::T
    function Segment(p1x, p1y, p2x, p2y)
        return new{Float64}(
            convert(Float64, p1x),
            convert(Float64, p1y),
            convert(Float64, p2x),
            convert(Float64, p2y),
        )
    end
end

function is_singular(s::Segment)
    if (s.p1x == s.p2x) || (s.p1y == s.p2y)
        return true
    else
        return false
    end
end

function intersect!(s1::Segment, s2::Segment, A::Matrix{Float64}, b::Vector{Float64})
    if is_singular(s1) || is_singular(s2)
        return false
    end
    A[1, 1] = s1.p2x - s1.p1x
    A[1, 2] = s2.p1x - s2.p2x
    A[2, 1] = s1.p2y - s1.p1y
    A[2, 2] = s2.p1y - s2.p2y
    b[1] = s2.p1x - s1.p1x
    b[2] = s2.p1y - s1.p1y
    sol = A \ b
    if (0 < sol[1] < 1) && (0 < sol[2] < 1)
        return true
    else
        return false
    end
end


function intersect(
    r1::Array{T},
    z1::Array{T},
    r2::Array{T},
    z2::Array{T},
) where {T<:AbstractFloat}
    A = zeros((2, 2))
    b = zeros(2)
    for i = 1:(length(r1) - 1)
        s1 = Segment(r1[i], z1[i], r1[i + 1], z1[i + 1])
        for j = 1:(length(r2) - 1)
            s2 = Segment(r2[j], z2[j], r2[j + 1], z2[j + 1])
            do_intersect = intersect!(s1, s2, A, b)
            if do_intersect
                return i
            end
        end
    end
    return length(r1)
end

intersect(integrator1::DenseIntegrator, integrator2::DenseIntegrator) =
    intersect(integrator1.r, integrator1.z, integrator2.r, integrator2.z)

function self_intersects(integrator::Sundials.IDAIntegrator)
    if length(integrator.p.data[:r]) < 2
        return false
    end
    r1 = integrator.p.data[:r][1:(end - 1)]
    z1 = integrator.p.data[:z][1:(end - 1)]
    r2 = integrator.p.data[:r][(end - 1):end]
    z2 = integrator.p.data[:z][(end - 1):end]
    index = intersect(r1, z1, r2, z2)
    if index < length(r1)
        return true
    else
        return false
    end
end

function reduce_integrators(
    integrators::Vector{<:Sundials.IDAIntegrator};
    n_timesteps = 1000,
    log = true,
)
    rs = Float64[]
    zs = Float64[]
    vrs = Float64[]
    vzs = Float64[]
    ns = Float64[]
    @info "Filtering intersections..."
    flush()
    #dense_integrators = DenseIntegrator.(integrators, n_timesteps = n_timesteps, log = log)
    dense_integrators = DenseIntegrator.(integrators, 0.1)#n_timesteps = n_timesteps, log = log)
    @showprogress for (i, integrator1) in enumerate(dense_integrators[1:(end - 1)])
        f(integ2) = intersect(integrator1, integ2)
        index = minimum(pmap(f, dense_integrators[(i + 1):end], batch_size=5))
        rs = vcat(rs, integrator1.r[1:index])
        zs = vcat(zs, integrator1.z[1:index])
        vrs = vcat(vrs, integrator1.vr[1:index])
        vzs = vcat(vzs, integrator1.vz[1:index])
        ns = vcat(ns, integrator1.n[1:index])
    end
    @info "Done"
    flush()
    return rs, zs, vrs, vzs, ns
end

function reduce_integrators_individual(integrators; n_timesteps = 1000, log = true)
    reduced_integrators = DenseIntegrator[]
    dense_integrators = DenseIntegrator.(integrators, n_timesteps = n_timesteps, log = log)
    for (i, integrator) in enumerate(dense_integrators[1:(end - 1)])
        f(integ2) = intersect(integrator, integ2)
        indices = pmap(f, dense_integrators[(i + 1):end])
        index = minimum(indices)
        reduced_integrator = DenseIntegrator(
            integrator.r[1:index],
            integrator.z[1:index],
            integrator.vr[1:index],
            integrator.vz[1:index],
            integrator.n[1:index],
        )
        push!(reduced_integrators, reduced_integrator)
    end
    return reduced_integrators
end
