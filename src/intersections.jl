import Base: intersect!, intersect
using LinearAlgebra, ProgressMeter
export Segment,
    intersect,
    reduce_integrators,
    self_intersects,
    interpolate_integrator,
    interpolate_integrators,
    get_intersection_times

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

struct DenseIntegrator{T<:Vector{<:AbstractFloat}}
    t::T
    r::T
    z::T
    vr::T
    vz::T
    n::T
end

intersect(integrator1::DenseIntegrator, integrator2::DenseIntegrator) =
    intersect(integrator1.r, integrator1.z, integrator2.r, integrator2.z)

function DenseIntegrator(integrator::Sundials.IDAIntegrator)
    return DenseIntegrator(
        integrator.sol.t[1:(end - 2)],
        integrator.p.data[:r],
        integrator.p.data[:z],
        integrator.p.data[:vr],
        integrator.p.data[:vz],
        integrator.p.data[:n],
    )
end

function self_intersects(integrator::Sundials.IDAIntegrator, r, z)
    if length(integrator.p.data[:r]) < 2
        return false
    end
    r1 = integrator.p.data[:r][1:(end - 1)]
    z1 = integrator.p.data[:z][1:(end - 1)]
    r2 = [integrator.p.data[:r][end], r] #integrator.p.data[:r][(end - 1):end]
    z2 = [integrator.p.data[:z][end], z] #integrator.p.data[:z][(end - 1):end]
    index = intersect(r1, z1, r2, z2)
    if index < length(r1)
        return true
    else
        return false
    end
end

function get_intersection_times(integrators::Vector{<:Sundials.IDAIntegrator})
    raw_integrators = DenseIntegrator.(integrators)
    intersection_times = Float64[]
    @showprogress for (i, raw_integrator) in enumerate(raw_integrators)
        if i == length(raw_integrators)
            # last one has absolute priority.
            index = length(raw_integrator.r)
        else
            f(integ2) = intersect(raw_integrator, integ2)
            index = minimum(pmap(f, raw_integrators[(i + 1):end], batch_size = 5))
        end
        push!(intersection_times, raw_integrator.t[index])
    end
    return intersection_times
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
        t_range = range(tmin, max_time, length = n_timesteps)
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
    line_width_dense = r_dense .* integrator.p.lwnorm
    density_dense =
        compute_density.(r_dense, z_dense, vr_dense, vz_dense, Ref(integrator.p))
    return DenseIntegrator(t_range, r_dense, z_dense, vr_dense, vz_dense, density_dense)
end

function interpolate_integrators(
    integrators::Vector{<:Sundials.IDAIntegrator};
    n_timesteps::Int = 1000,
    max_times = nothing,
    log = true,
)
    ret = DenseIntegrator[]
    for (integrator, max_time) in zip(integrators, max_times)
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
    return ret
end

function reduce_integrators(
    integrators::Vector{<:Sundials.IDAIntegrator};
    max_times = nothing,
)
    time_indices =
        searchsorted_nearest.([integrator.sol.t for integrator in integrators], max_times)
    r = reduce(
        vcat,
        [
            integrator.p.data[:r][1:index]
            for (index, integrator) in zip(time_indices, integrators)
        ],
    )
    z = reduce(
        vcat,
        [
            integrator.p.data[:z][1:index]
            for (index, integrator) in zip(time_indices, integrators)
        ],
    )
    vr = reduce(
        vcat,
        [
            integrator.p.data[:vr][1:index]
            for (index, integrator) in zip(time_indices, integrators)
        ],
    )
    vz = reduce(
        vcat,
        [
            integrator.p.data[:vz][1:index]
            for (index, integrator) in zip(time_indices, integrators)
        ],
    )
    n = reduce(
        vcat,
        [
            integrator.p.data[:n][1:index]
            for (index, integrator) in zip(time_indices, integrators)
        ],
    )
    return r, z, vr, vz, n
end

function reduce_integrators(integrators::Vector{<:DenseIntegrator})
    r = reduce(vcat, [integrator.r for integrator in integrators])
    z = reduce(vcat, [integrator.z for integrator in integrators])
    vr = reduce(vcat, [integrator.vr for integrator in integrators])
    vz = reduce(vcat, [integrator.vz for integrator in integrators])
    n = reduce(vcat, [integrator.n for integrator in integrators])
    return r, z, vr, vz, n
end

function reduce_integrators(integrators::Vector{<:DenseIntegrator})
    r = reduce(vcat, [integrator.r for integrator in integrators])
    z = reduce(vcat, [integrator.z for integrator in integrators])
    vr = reduce(vcat, [integrator.vr for integrator in integrators])
    vz = reduce(vcat, [integrator.vz for integrator in integrators])
    n = reduce(vcat, [integrator.n for integrator in integrators])
    return r, z, vr, vz, n
end
