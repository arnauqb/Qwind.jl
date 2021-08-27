using Qwind, Test

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

@testset "Test interpolation grid" begin

    r_range = Float64[]
    z_range = Float64[]
    for r in range(10,100, length=50)
        for z in range(10,100, length=50)
            push!(r_range, r)
            push!(z_range, z)
        end
    end
    hull = Hull(
        [0.01, 200, 0.01, 200],
        [0.01, 200, 200, 0.01],
    );
    r0_range = collect(range(0.1, 100.0, step = 10));
    func(r, z) = 1e12 * exp(-z / 100) / r;
    n_range = func.(r_range, z_range);
    grid = DensityGrid(
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

