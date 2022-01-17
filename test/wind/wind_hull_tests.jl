using Qwind, Test

@testset "Construct wind hull" begin
    n = 10000
    # make a triangle
    r_randoms = 100 .* rand(n)
    z_randoms = 100 .* rand(n)
    rs = Float64[]
    zs = Float64[]
    for (rp, zp) in zip(r_randoms, z_randoms)
        if rp >= zp
            push!(rs, rp)
            push!(zs, zp)
        end
    end
    @test length(rs) == length(zs)
    hull = Hull(rs, zs)
    @test Qwind.is_point_in_wind(hull, [0, 10]) == false
    @test Qwind.is_point_in_wind(hull, [50, 10]) == true
    @test Qwind.is_point_in_wind(hull, [90, 80]) == true
    @test Qwind.is_point_in_wind(hull, [40, 30]) == true
    @test Qwind.is_point_in_wind(hull, [200, 0]) == false
end
