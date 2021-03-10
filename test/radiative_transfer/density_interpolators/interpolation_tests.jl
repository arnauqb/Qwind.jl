using Qwind, Test

@testset "Reduce line" begin
    r = [1, 2, 3, 4, 5]
    z = [1, 2, 3, 4, 3]
    n = [10, 20, 30, 40, 50]
    rr, zz, nn = Qwind.reduce_line(r, z, n)
    @test rr == r[1:(end - 1)]
    @test zz == z[1:(end - 1)]
    @test nn == n[1:(end - 1)]
end

@testset "Get  spatial grid" begin
    r = [1.0, 2.0, 30.0, 4.0, 105.0]
    z = [1.0, 2.0, 3.0, 4.0, 3.0]
    r0s = [1.0, 2.0, 3.0]
    auto_r, auto_z = Qwind.get_spatial_grid(r, z, r0s, "auto", 20)
    @test auto_r[1:3] == r0s
    for i = 1:(length(auto_r) - 1)
        @test auto_r[i + 1] > auto_r[i]
    end
    @test auto_r[end] == 103.0
    @test auto_z[1] == 1.0
    @test auto_z[end] == 4
    @test length(auto_z) == 19 # because one is left for z = 0
    for i = 1:(length(auto_z) - 1)
        @test auto_z[i + 1] > auto_z[i]
    end
end

#@testset "Test empty interpolation grid" begin
#    r_range = range(0, 1000, length = 100)
#    z_range = range(0, 1000, length = 100)
#    empty_grid = Qwind.construct_interpolation_grid(10, 20)
#    @test empty_grid.nr == 10
#    @test empty_grid.nz == 20
#    @test empty_grid.r_range == [0.0, 0.0]
#    @test empty_grid.z_range == [0.0, 0.0]
#    for (i, r) in enumerate(r_range)
#        for (j, z) in enumerate(z_range)
#            @test empty_grid.interpolator(r, z) == 1e2
#        end
#    end
#end

@testset "Construct wind hull" begin
    rr = collect(range(1.0, 100, length = 10))
    rr = vcat(rr, 100.0 .* ones(10))
    zz = collect(range(1.0, 100, length = 10))
    zz = vcat(zz, range(1.0, 100, length=10))
    #rr = vcat(rr, range(0, 100, length = 20))
    #zz = vcat(zz, zeros(20))
    @test length(rr) == length(zz)
    r0s = collect(range(1.0, 100.0, length = 5))
    hull = Qwind.construct_wind_hull(rr, zz, r0s)
    @test Qwind.is_point_in_wind(hull, [0, 10]) == false
    @test Qwind.is_point_in_wind(hull, [50, 0]) == false
    @test Qwind.is_point_in_wind(hull, [99, 99]) == true
    @test Qwind.is_point_in_wind(hull, [99, 100]) == false
    @test Qwind.is_point_in_wind(hull, [99, 98]) == true
    @test Qwind.is_point_in_wind(hull, [40, 30]) == true
    @test Qwind.is_point_in_wind(hull, [200, 0]) == false
end

@testset "Test interpolation grid" begin

    r_range = Float64[]
    z_range = Float64[]
    for r in range(10,100, length=50)
        for z in range(10,100, length=50)
            push!(r_range, r)
            push!(z_range, z)
        end
    end
    hull = Qwind.construct_wind_hull(
        [0.01, 200, 0.01, 200],
        [0.01, 200, 200, 0.01],
        [0.01, 200.0, 0.01, 200.0]
    );
    r0_range = collect(range(0.1, 100.0, step = 10));
    func(r, z) = 1e12 * exp(-z / 100) / r;
    n_range = func.(r_range, z_range);
    grid = Qwind.construct_interpolation_grid(
        r_range,
        z_range,
        n_range,
        r0_range,
        hull,
        nr = 50,
        nz = 50,
    );

    @testset "Test gridpoints" begin
        # Grid points should be exact?
        for (i, r) in enumerate(grid.r_range)
            for (j, z) in enumerate(grid.z_range)
                @test grid.grid[i,j] ≈ func(r,z) rtol = 0.1
            end
        end
    end

    @testset "Test interpolation" begin
        for r in range(grid.r_range[1], grid.r_range[end], length=100)
            for z in range(grid.z_range[1], grid.z_range[end], length=100)
                @test func(r,z) ≈ grid.interpolator(r,z) rtol = 0.1
            end
        end
    end

end


