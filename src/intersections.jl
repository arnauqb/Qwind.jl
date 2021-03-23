import Base: intersect!, intersect
using LinearAlgebra, ProgressMeter
export Segment,
    intersect,
    reduce_integrators,
    self_intersects,
    interpolate_integrator,
    interpolate_integrators,
    get_intersection_times,
    load_trajectories

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
    n1::Array{T},
    v1::Array{T},
    r2::Array{T},
    z2::Array{T},
    n2::Array{T},
    v2::Array{T},
) where {T<:AbstractFloat}
    A = zeros((2, 2))
    b = zeros(2)
    for i = 1:(length(r1) - 1)
        s1 = Segment(r1[i], z1[i], r1[i + 1], z1[i + 1])
        m1 = n1[i] * v1[i]^2
        for j = 1:(length(r2) - 1)
            s2 = Segment(r2[j], z2[j], r2[j + 1], z2[j + 1])
            m2 = n2[j] * v2[j]^2
            do_intersect = intersect!(s1, s2, A, b)
            if do_intersect
                if m2 > m1 # terminate when intersecting has higher momentum
                    return i
                else
                    continue
                end
            end
        end
    end
    return length(r1)
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

function intersect(integrator1::DenseIntegrator, integrator2::DenseIntegrator)
    v1 = sqrt.(integrator1.vr .^ 2 + integrator1.vz .^ 2)
    v2 = sqrt.(integrator2.vr .^ 2 + integrator2.vz .^ 2)
    intersect(
        integrator1.r,
        integrator1.z,
        integrator1.n,
        v1,
        integrator2.r,
        integrator2.z,
        integrator2.n,
        v2,
    )
end

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

function load_trajectory(tdata::Dict)
    return DenseIntegrator(tdata["t"], tdata["r"], tdata["z"], tdata["vr"], tdata["vz"], tdata["n"])
end

function load_trajectories(tsdata::Dict)
    t_ids = sort(parse.(Int, keys(tsdata)))
    t_ids = string.(t_ids)
    ret = DenseIntegrator[]
    for i in t_ids
        tdata = tsdata[i]
        trajectory = load_trajectory(tdata)
        push!(ret, trajectory)
    end
    return ret
end

function load_trajectories(h5_path::String, it_num)
    it_name = @sprintf "iteration_%03d" it_num
    tsdata = h5open(h5_path, "r") do file
        read(file, it_name * "/trajectories")
    end
    return load_trajectories(tsdata)
end

function load_trajectories(h5_path::String)
    it_keys = h5open(h5_path, "r") do file
        keys(read(file))
    end
    it_nums = [parse(Int, split(key, "_")[end]) for key in it_keys]
    return load_trajectories(h5_path, maximum(it_nums))
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

function slice_integrator(integrator::DenseIntegrator; in = 1, fi = nothing)
    (fi === nothing) && (fi = length(integrator.t))
    return DenseIntegrator(
        integrator.t[in:fi],
        integrator.r[in:fi],
        integrator.z[in:fi],
        integrator.vr[in:fi],
        integrator.vz[in:fi],
        integrator.n[in:fi],
    )
end

function get_intersection_times(integrators::Vector{<:DenseIntegrator})
    intersection_indices = [length(integrator.t) for integrator in integrators]
    @info "Calculating trajectory intersections..."
    n_procs_to_use = min(10, nprocs() - 1)
    worker_pool = WorkerPool(collect(1:n_procs_to_use))
    @showprogress for (i, integrator) in enumerate(integrators)
        integrators_sliced = [
            slice_integrator(integrator, fi = index)
            for (integrator, index) in zip(integrators, intersection_indices)
        ]
        index = length(integrator.t)
        integrator_indcs = [j for j in 1:length(integrators) if j != i]
        #@sync @distributed min for integrator2 = integrators_sliced[integrator_indcs]
        #    intersect(integrator, integrator2)
        #end
        f(j) = intersect(integrator, integrators_sliced[j])
        #integrator_indcs = [j for j in 1:length(integrators) if j != i]
        #batch = Int(floor(length(integrator_indcs) // n_procs_to_use))
        index = minimum(pmap(f, worker_pool, integrator_indcs, batch_size = 50))
        #println(worker_pool)
        #index = minimum(pmap(f, integrator_indcs, batch_size = 50))
        intersection_indices[i] = index
    end
    intersection_times = [
        integrator.t[index]
        for (integrator, index) in zip(integrators, intersection_indices)
    ]
    return intersection_times
end

get_intersection_times(integrators::Vector{<:Sundials.IDAIntegrator}) =
    get_intersection_times(DenseIntegrator.(integrators))

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

