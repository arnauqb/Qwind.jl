using LinearAlgebra, ProgressMeter, DataStructures
export Segment,
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

min_x(segment::Segment) = min(segment.p1x, segment.p2x)
min_y(segment::Segment) = min(segment.p1y, segment.p2y)
max_x(segment::Segment) = max(segment.p1x, segment.p2x)
max_y(segment::Segment) = max(segment.p1y, segment.p2y)

function is_singular(s::Segment)
    if (s.p1x == s.p2x) || (s.p1y == s.p2y)
        return true
    else
        return false
    end
end

struct Intersection{S<:String, I<:Int, T<:Float64}
    type::S
    id1::I
    id2::I
    t1::T
    t2::T
end

function Base.isless(i1::Intersection{S, I, T}, i2::Intersection{S, I, T}) where {S<:String, I<:Int, T<:Float64}
    if (i1.t1 < i2.t1)
        return true
    elseif i1.t1 == i2.t1
        if i1.t2 < i2.t2
            return true
        elseif i1.t2 == i2.t2
            return i1.id1 < i2.id1
        end
    end
    return false
end

function trivial_miss(s1::Segment, s2::Segment)
    if max_x(s1) < min_x(s2)
        return true
    elseif max_x(s2) < min_x(s1)
        return true
    elseif max_y(s1) < min_y(s2)
        return true
    elseif max_y(s2) < min_y(s1)
        return true
    end
    return false
end

function Base.intersect!(s1::Segment, s2::Segment, A::Matrix{Float64}, b::Vector{Float64})
    if is_singular(s1) || is_singular(s2)
        return false
    end
    if trivial_miss(s1, s2)
        return false
    end
    if s1 == s2
        return true
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


function Base.intersect!(
    intersections_dict::Dict{Int,SortedSet{Intersection}},
    id1::Int,
    t1::Array{T},
    r1::Array{T},
    z1::Array{T},
    n1::Array{T},
    v1::Array{T},
    id2::Int,
    t2::Array{T},
    r2::Array{T},
    z2::Array{T},
    n2::Array{T},
    v2::Array{T},
) where {T<:AbstractFloat}
    A = zeros((2, 2))
    b = zeros(2)
    for i = 1:(length(r1) - 1)
        s1 = Segment(r1[i], z1[i], r1[i + 1], z1[i + 1])
        m1 = sign(v1[i]) * n1[i] * v1[i]^2
        for j = 1:(length(r2) - 1)
            s2 = Segment(r2[j], z2[j], r2[j + 1], z2[j + 1])
            m2 = sign(v2[j]) * n2[j] * v2[j]^2
            do_intersect = intersect!(s1, s2, A, b)
            if do_intersect
                if m1 > m2
                    intersection = Intersection("winner", id1, id2, t1[i], t2[j])
                    insert_intersection!(intersections_dict, id1, intersection)
                    intersection = Intersection("loser", id2, id1, t2[j], t1[i])
                    insert_intersection!(intersections_dict, id2, intersection)
                else
                    intersection = Intersection("winner", id2, id1, t2[j], t1[i])
                    insert_intersection!(intersections_dict, id2, intersection)
                    intersection = Intersection("loser", id1, id2, t1[i], t2[j])
                    insert_intersection!(intersections_dict, id1, intersection)
                end
            end
        end
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
    id::Int
    t::T
    r::T
    z::T
    vr::T
    vz::T
    n::T
end

function Base.intersect!(
    intersections_dict::Dict{Int,SortedSet{Intersection}},
    integrator1::DenseIntegrator,
    integrator2::DenseIntegrator,
)
    v1 = sign.(integrator1.vz) .* sqrt.(integrator1.vr .^ 2 + integrator1.vz .^ 2)
    v2 = sign.(integrator2.vz) .* sqrt.(integrator2.vr .^ 2 + integrator2.vz .^ 2)
    intersect!(
        intersections_dict,
        integrator1.id,
        integrator1.t,
        integrator1.r,
        integrator1.z,
        integrator1.n,
        v1,
        integrator2.id,
        integrator2.t,
        integrator2.r,
        integrator2.z,
        integrator2.n,
        v2,
    )
end

function DenseIntegrator(integrator::Sundials.IDAIntegrator)
    return DenseIntegrator(
        integrator.p.id,
        integrator.sol.t[1:(end - 2)],
        integrator.p.data[:r],
        integrator.p.data[:z],
        integrator.p.data[:vr],
        integrator.p.data[:vz],
        integrator.p.data[:n],
    )
end

function load_trajectory(tdata::Dict)
    return DenseIntegrator(
        tdata["id"],
        tdata["t"],
        tdata["r"],
        tdata["z"],
        tdata["vr"],
        tdata["vz"],
        tdata["n"],
    )
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
    r2 = [integrator.p.data[:r][end], r] 
    z2 = [integrator.p.data[:z][end], z] 
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
        integrator.id,
        integrator.t[in:fi],
        integrator.r[in:fi],
        integrator.z[in:fi],
        integrator.vr[in:fi],
        integrator.vz[in:fi],
        integrator.n[in:fi],
    )
end

function get_intersections(integrators::Vector{<:DenseIntegrator})
    intersections_dict = Dict{Int,SortedSet{Intersection}}()
    for i = 1:length(integrators)
        intersections = SortedSet{Intersection}()
        for j = (i+1):length(integrators)
            intersect!(intersections_dict, integrators[i], integrators[j])
        end
    end
    # Turn Sets to vectors
    ret = Dict{Int, Vector{Intersection{String, Int, Float64}}}()
    for (key, intersections) in intersections_dict
        ret[key] = collect(intersections)
    end
    return ret
end

function insert_intersection!(dict::Dict{Int, SortedSet{Intersection}}, key, intersection)
    if !haskey(dict, key)
        dict[key] = SortedSet([intersection])
    else
        push!(dict[key], intersection)
    end
end

function insert_intersection!(dict, key, intersection)
    if !haskey(dict, key)
        dict[key] = [intersection]
    else
        push!(dict[key], intersection)
    end
end

function delete_intersection!(
    dict::Dict{Int,Vector{Intersection}},
    key::Int,
    intersection::Intersection,
)
    delete!(dict[key], intersection)
    if length(dict[key]) == 0
        delete!(dict, key)
    end
end

function clear_remaining_trajectory_intersections!(intersections_dict, traj_id, max_time)
    intersections = intersections_dict[traj_id]
    todelete_dict = DefaultDict{Int, Set{Int}}(() -> Set{Int}())
    for (i, intersection) in enumerate(intersections)
        if intersection.t1 >= max_time
            push!(todelete_dict[traj_id], i)
            for (j, intersection2) in enumerate(intersections_dict[intersection.id2])
                if intersection2.t2 == intersection.t1
                    push!(todelete_dict[intersection.id2], j)
                end
            end
        end
    end
    for (key, todelete) in todelete_dict
        deleteat!(intersections_dict[key], sort(collect(todelete)))
        if length(intersections_dict[key]) == 0
            delete!(intersections_dict, key)
        end
    end
end

function resolve_intersections!(intersection_times, intersections_dict)
    n_inters = length(intersections_dict)
    while length(intersections_dict) > 0
        all_keys = keys(intersections_dict)
        stuck = true
        for traj_id in all_keys
            intersections = intersections_dict[traj_id]
            found = false
            for intersection in intersections
                if intersection.type == "loser"
                    break
                end
                stuck = false
                loser_id = intersection.id2
                max_time = intersection.t2
                intersection_times[loser_id] = min(intersection_times[loser_id], max_time)
                clear_remaining_trajectory_intersections!(
                    intersections_dict,
                    loser_id,
                    max_time,
                )
                found = true
                break# to clear index
            end
            (found == true) && break
        end
        if stuck
            # this means that they all have a "losing intersection" as their first.
            # Execute the innermost one
            key = minimum(all_keys)
            intersection = first(intersections_dict[key])
            loser_id = intersection.id2
            max_time = intersection.t2
            intersection_times[loser_id] = min(intersection_times[loser_id], max_time)
            clear_remaining_trajectory_intersections!(
                intersections_dict,
                loser_id,
                max_time,
            )
        end
    end
    return intersection_times
end

function get_intersection_times(integrators::Vector{<:DenseIntegrator})
    intersection_times =
        Dict(integrator.id => integrator.t[end] for integrator in integrators)
    intersections_dict = get_intersections(integrators)
    resolve_intersections!(intersection_times, intersections_dict)
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
    return DenseIntegrator(
        integrator.p.id,
        t_range,
        r_dense,
        z_dense,
        vr_dense,
        vz_dense,
        density_dense,
    )
end

function interpolate_integrators(
    integrators::Vector{<:Sundials.IDAIntegrator};
    n_timesteps::Int = 1000,
    max_times = nothing,
    log = true,
)
    ret = DenseIntegrator[]
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
    return ret
end

function reduce_integrators(
    integrators::Vector{<:Sundials.IDAIntegrator};
    max_times = nothing,
)
    max_times = [max_times[integrator.id] for integrator in integrators]
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

