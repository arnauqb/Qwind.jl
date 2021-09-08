using Qwind, Test

@testset "save streamlines properties" begin

    # Think about solving an entire system where all the mass is lost?
    t = collect(range(1, 10, length=100))
    r = collect(range(10, 1000, length=100))
    z = collect(range(20, 2000, length=100))
    vr = collect(range(10, 1000, length=100))
    vphi = collect(range(10, 1000, length=100))
    vz = collect(range(20, 2000, length=100))
    n = collect(range(1e8, 1e9, length=100));

    integrators = []
    for i in 1:100
        push!(integrators, Qwind.DenseIntegrator(i, t, r, z, vr, vphi, vz, n))
    end

    bh = BlackHole(1e8 * M_SUN, 0.5, 0.0)

    wind_properties = Qwind.create_wind_properties(integrators, bh);
    vt = sqrt.(vr .^2 + vz.^2);
    mass_loss = π * (r[end]^2 - r[1]^2) * n[1] * M_P / 0.61 * vt[1] * C;

    @test wind_properties["mass_loss"] ≈ mass_loss




end

