import Base: intersect
using LinearAlgebra
export Segment, intersect

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

function intersect(s1::Segment, s2::Segment)
    A = [s1.p2x-s1.p1x s2.p1x-s2.p2x; s1.p2y-s1.p1y s2.p1y-s2.p2y]
    b = [s2.p1x - s1.p1x, s2.p1y - s1.p1y]
    sol = A \ b
    if (0 < sol[1] < 1) && (0 < sol[2] < 1)
        sol = [s1.p1x + sol[1] * (s1.p2x - s1.p1x), s1.p1y + sol[1] * (s1.p2y - s1.p1y)]
        return true, sol
    else
        return false, [Inf, Inf]
    end
end


function intersect(
    r1::Array{T},
    z1::Array{T},
    r2::Array{T},
    z2::Array{T},
) where {T<:AbstractFloat}
    for i = 1:length(r1) - 1
        s1 = Segment(r1[i], z1[i], r1[i + 1], z1[i + 1])
        for j = 1:length(r2)-1
            s2 = Segment(r2[j], z2[j], r2[j + 1], z2[j + 1])
            do_intersect, intersection = intersect(s1, s2)
            if do_intersect
                return do_intersect, intersection
            end
        end
    end
    return false, [Inf, Inf]
end
