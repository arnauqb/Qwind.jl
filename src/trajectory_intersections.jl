using LinearAlgebra, ProgressMeter, DataStructures
export Segment, Intersection, reduce_trajectories, get_intersection_times, load_trajectories

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

struct Intersection{S<:String,I<:Int,T<:Float64}
    type::S
    id1::I
    id2::I
    t1::T
    t2::T
end

function Base.isless(
    i1::Intersection{S,I,T},
    i2::Intersection{S,I,T},
) where {S<:String,I<:Int,T<:Float64}
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
    try
        sol = A \ b
        if (0 < sol[1] < 1) && (0 < sol[2] < 1)
            return true
        else
            return false
        end
    catch e
        if isa(e, LinearAlgebra.SingularException)
            # this happens if they are parallel
            return false
        else
            throw(DomainError)
        end
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
    escaped1::Bool,
    id2::Int,
    t2::Array{T},
    r2::Array{T},
    z2::Array{T},
    n2::Array{T},
    v2::Array{T},
    escaped2::Bool,
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
                if escaped1
                    if escaped2
                        if m1 > m2
                            # 1 wins
                            intersection = Intersection("winner", id1, id2, t1[i], t2[j])
                            insert_intersection!(intersections_dict, id1, intersection)
                            intersection = Intersection("loser", id2, id1, t2[j], t1[i])
                            insert_intersection!(intersections_dict, id2, intersection)
                        else
                            # 2 wins
                            intersection = Intersection("winner", id2, id1, t2[j], t1[i])
                            insert_intersection!(intersections_dict, id2, intersection)
                            intersection = Intersection("loser", id1, id2, t1[i], t2[j])
                            insert_intersection!(intersections_dict, id1, intersection)
                        end
                    else
                        # 1 wins
                        intersection = Intersection("winner", id1, id2, t1[i], t2[j])
                        insert_intersection!(intersections_dict, id1, intersection)
                        intersection = Intersection("loser", id2, id1, t2[j], t1[i])
                        insert_intersection!(intersections_dict, id2, intersection)
                    end
                else
                    if escaped2
                        # 2 wins
                        intersection = Intersection("winner", id2, id1, t2[j], t1[i])
                        insert_intersection!(intersections_dict, id2, intersection)
                        intersection = Intersection("loser", id1, id2, t1[i], t2[j])
                        insert_intersection!(intersections_dict, id1, intersection)
                    else
                        if m1 > m2
                            # 1 wins
                            intersection = Intersection("winner", id1, id2, t1[i], t2[j])
                            insert_intersection!(intersections_dict, id1, intersection)
                            intersection = Intersection("loser", id2, id1, t2[j], t1[i])
                            insert_intersection!(intersections_dict, id2, intersection)
                        else
                            # 2 wins
                            intersection = Intersection("winner", id2, id1, t2[j], t1[i])
                            insert_intersection!(intersections_dict, id2, intersection)
                            intersection = Intersection("loser", id1, id2, t1[i], t2[j])
                            insert_intersection!(intersections_dict, id1, intersection)
                        end
                    end
                end
            end
        end
    end
end

function Base.intersect(
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

struct Trajectory{T<:Vector{<:AbstractFloat}}
    id::Int
    t::T
    r::T
    z::T
    vr::T
    vphi::T
    vz::T
    n::T
    escaped::Bool
end

function Base.intersect!(
    intersections_dict::Dict{Int,SortedSet{Intersection}},
    trajectory1::Trajectory,
    trajectory2::Trajectory,
)
    v1 = @. sign(trajectory1.vz) * sqrt.(trajectory1.vr^2 + trajectory1.vz^2)
    v2 = @. sign(trajectory2.vz) * sqrt.(trajectory2.vr^2 + trajectory2.vz^2)
    n_intersections = intersect!(
        intersections_dict,
        trajectory1.id,
        trajectory1.t,
        trajectory1.r,
        trajectory1.z,
        trajectory1.n,
        v1,
        trajectory1.escaped,
        trajectory2.id,
        trajectory2.t,
        trajectory2.r,
        trajectory2.z,
        trajectory2.n,
        v2,
        trajectory2.escaped,
    )
    return n_intersections
end

function Trajectory(integrator::Sundials.IDAIntegrator)
    return Trajectory(
        integrator.p.id,
        integrator.sol.t[1:(end - 2)],
        integrator.p.data[:r],
        integrator.p.data[:z],
        integrator.p.data[:vr],
        integrator.p.data[:vphi],
        integrator.p.data[:vz],
        integrator.p.data[:n],
        escaped(integrator),
    )
end

function load_trajectory(tdata::Dict)
    return Trajectory(
        tdata["id"],
        tdata["t"],
        tdata["r"],
        tdata["z"],
        tdata["vr"],
        tdata["vphi"],
        tdata["vz"],
        tdata["n"],
        tdata["escaped"],
    )
end

function load_trajectories(tsdata::Dict)
    t_ids = sort(parse.(Int, keys(tsdata)))
    t_ids = string.(t_ids)
    ret = Trajectory[]
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

function slice_trajectory(trajectory::Trajectory; in = 1, fi = nothing)
    (fi === nothing) && (fi = length(trajectory.t))
    return Trajectory(
        trajectory.id,
        trajectory.t[in:fi],
        trajectory.r[in:fi],
        trajectory.z[in:fi],
        trajectory.vr[in:fi],
        trajectory.vphi[in:fi],
        trajectory.vz[in:fi],
        trajectory.n[in:fi],
        trajectory.escaped,
    )
end

function get_intersections(trajectories::Vector{<:Trajectory})
    intersections_dict = Dict{Int,SortedSet{Intersection}}()
    n_inters = 0
    for i = 1:length(trajectories)
        intersections = SortedSet{Intersection}()
        for j = (i + 1):length(trajectories)
            intersect!(intersections_dict, trajectories[i], trajectories[j])
        end
    end
    # Turn Sets to vectors
    ret = Dict{Int,Vector{Intersection{String,Int,Float64}}}()
    for (key, intersections) in intersections_dict
        ret[key] = collect(intersections)
    end
    return ret
end

function insert_intersection!(dict::Dict{Int,SortedSet{Intersection}}, key, intersection)
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
    todelete_dict = DefaultDict{Int,Set{Int}}(() -> Set{Int}())
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

function resolve_intersections(intersection_times, intersections_dict)
    n_inters = length(intersections_dict)
    merging_streamlines = Dict()
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
                winner_id = intersection.id1
                max_time = intersection.t2
                winner_time = intersection.t1
                if max_time < intersection_times[loser_id]
                    intersection_times[loser_id] = max_time
                    merging_streamlines[loser_id] = Dict(
                        "winner_id" => winner_id,
                        "winner_time" => winner_time,
                        "loser_time" => max_time,
                    )
                end
                clear_remaining_trajectory_intersections!(
                    intersections_dict,
                    loser_id,
                    max_time,
                )
                found = true
                break # to clear index
            end
            (found == true) && break
        end
        if stuck
            # this means that they all have a "losing intersection" as their first.
            # Execute the innermost one
            key = minimum(all_keys)
            intersection = first(intersections_dict[key])
            loser_id = intersection.id2
            winner_id = intersection.id1
            max_time = intersection.t2
            if max_time < intersection_times[loser_id]
                intersection_times[loser_id] = max_time
                merging_streamlines[loser_id] = winner_id
            end
            clear_remaining_trajectory_intersections!(
                intersections_dict,
                loser_id,
                max_time,
            )
        end
    end
    return intersection_times, merging_streamlines
end

function get_intersection_times(integrators::Vector{<:Trajectory})
    intersection_times =
        Dict(integrator.id => integrator.t[end] for integrator in integrators)
    intersections_dict = get_intersections(integrators)
    intersection_times, merging_streamlines =
        resolve_intersections(intersection_times, intersections_dict)
    return intersection_times, merging_streamlines
end

get_intersection_times(integrators::Vector{<:Sundials.IDAIntegrator}) =
    get_intersection_times(Trajectory.(integrators))


# Line widths calculations

function get_distance_between_segments(
    s1::Segment,
    s2::Segment,
    A::Matrix{Float64},
    b::Vector{Float64},
)
    s1_mid_x = s1.p1x + (s1.p2x - s1.p1x) / 2.0
    s1_mid_y = s1.p1y + (s1.p2y - s1.p1y) / 2.0
    s1_vec_x = s1.p2x - s1.p1x
    s1_vec_y = s1.p2y - s1.p1y
    s2_vec_x = s2.p2x - s2.p1x
    s2_vec_y = s2.p2y - s2.p1y
    s1_perp_vec_x = 1.0
    s1_perp_vec_y = -s1_vec_x / s1_vec_y
    A[1, 1] = s1_perp_vec_x
    A[1, 2] = -s2_vec_x
    A[2, 1] = s1_perp_vec_y
    A[2, 2] = -s2_vec_y
    b[1] = s2.p1x - s1_mid_x
    b[2] = s2.p1y - s1_mid_y
    sol = A \ b
    if (0 < sol[2] < 1)
        int_x = s2.p1x + sol[2] * s2_vec_x
        int_y = s2.p1y + sol[2] * s2_vec_y
        distance = sqrt((int_x - s1_mid_x)^2 + (int_y - s1_mid_y)^2)
    else
        distance = Inf
    end
    return distance
end

function turn_streamline_to_segments(sl)
    ret = Array{Segment{typeof(sl.r[1])}}(undef, length(sl.r) - 1)
    for i = 1:(length(sl.r) - 1)
        seg = Segment(sl.r[i], sl.z[i], sl.r[i + 1], sl.z[i + 1])
        ret[i] = seg
    end
    return ret
end


function get_width_pairs(centre, left, right, A, b)
    centre_segs = turn_streamline_to_segments(centre)
    left_segs = turn_streamline_to_segments(left)
    right_segs = turn_streamline_to_segments(right)
    ret = []
    for centre_seg in centre_segs
        lws = [Inf, Inf]
        for left_seg in left_segs
            d = get_distance_between_segments(centre_seg, left_seg, A, b)
            if d < lws[1]
                lws[1] = d
            end
        end
        for right_seg in right_segs
            d = get_distance_between_segments(centre_seg, right_seg, A, b)
            if d < lws[2]
                lws[2] = d
            end
        end
        push!(ret, lws)
    end
    return ret
end
