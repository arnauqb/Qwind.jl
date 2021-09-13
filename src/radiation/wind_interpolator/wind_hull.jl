import ConcaveHull
import ConcaveHull: Hull
import Sundials
export Hull

function Hull(r::Vector{Float64}, z::Vector{Float64})
    points = [[r[i], z[i]] for i = 1:length(r)]
    hull = ConcaveHull.concave_hull(points)
    return hull
end


function Hull(streamlines)
    r = reduce(vcat, [streamline.r[1:10:end] for streamline in streamlines])
    z = reduce(vcat, [streamline.z[1:10:end] for streamline in streamlines])
    @info "Constructing wind hull"
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
