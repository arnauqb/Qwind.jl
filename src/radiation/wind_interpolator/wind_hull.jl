import ConcaveHull
import ConcaveHull: Hull
import Sundials
export Hull

function Hull(r::Vector{Float64}, z::Vector{Float64})
    points = [[r[i], z[i]] for i = 1:length(r)]
    hull = ConcaveHull.concave_hull(points)
    return hull
end

function reduce_streamlines(streamlines; rtol = 0.1, atol=1e-4)
    r = Float64[]
    z = Float64[]
    for streamline in streamlines
        current_d = 0.0
        push!(r, streamline.r[1])
        push!(z, streamline.z[1])
        for i in 2:length(streamline.r)
            new_d = sqrt((streamline.r[i] - streamline.r[1])^2 + (streamline.z[i] - streamline.z[1])^2)
            rtol_condition = (abs((new_d - current_d) / new_d)) > rtol
            atol_condition = abs(new_d - current_d) > atol
            if rtol_condition && atol_condition
                push!(r, streamline.r[i])
                push!(z, streamline.z[i])
                current_d = new_d
            end
        end
    end
    # add r0 line
    r0s = range(streamlines[1].r[1], streamlines[length(streamlines)].r[1], length=2000);
    z0s = range(streamlines[1].z[1], streamlines[1].z[1], length=2000);
    r = vcat(r, r0s)
    z = vcat(z, z0s)
    return r, z
end

function Hull(streamlines; rtol=1e-3, atol=1e-2)
    @info "Constructing wind hull"
    r, z = reduce_streamlines(streamlines, rtol=rtol, atol=atol)
    hull = Hull(r, z)
    if !hull.converged
        error("Hull did not converge!")
    end
    @info "Done"
    return hull
end

function Hull(hull_data::Dict)
    vs_r = hull_data["vertices_r"]
    vs_z = hull_data["vertices_z"]
    vertices = [[r, z] for (r, z) in zip(vs_r, vs_z)]
    Hull(vertices, hull_data["k"], hull_data["converged"])
end


function Hull(h5_path::String, it_num)
    it_name = @sprintf "iteration_%03d" it_num
    grid_data = h5open(h5_path, "r") do file
        read(file, it_name * "/wind_hull")
    end
    return Hull(grid_data)
end

function Hull(h5_path::String)
    it_keys = h5open(h5_path, "r") do file
        keys(read(file))
    end
    it_nums = [parse(Int, split(key, "_")[end]) for key in it_keys]
    return Hull(h5_path, maximum(it_nums))
end

Hull() = ConcaveHull.concave_hull([[1e8, 1e8], [-1e8, -1e8], [1e8, -1e8], [-1e8, 1e8]])

function is_point_in_wind(hull::ConcaveHull.Hull, point)
    return ConcaveHull.in_hull(point, hull)
end
is_point_in_wind(hull::ConcaveHull.Hull, r, z) = is_point_in_wind(hull, [r, z])
