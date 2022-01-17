using Qwind, Test, HDF5


@testset "Save streamlines properties" begin

    # Think about solving an entire system where all the mass is lost?
    t = [1.0, 2.0]
    r = [100.0, 200.0]
    z = [0.0, 1.0]
    vr1 = [0.1, 0.3]
    vz1 = [0.1, 0.2]
    vr2 = [0.005, 0.007]
    vz2 = [0.001, 0.0015]
    vphi = [0.1, 0.1]
    n = [1e8, 1e9]
    widths = [1.0, 2.0]
    streamlines = []
    push!(streamlines, Streamline(1, t, r, z, vr1, vphi, vz1, n, widths, true)) # this escapes
    push!(streamlines, Streamline(2, t, r, z, vr2, vphi, vz2, n, widths, false)) # this doesn't
    streamlines = Streamlines(streamlines)
    bh = BlackHole(1e8 * M_SUN, 0.5, 0.0)
    @test escaped(streamlines[1]) == true
    @test escaped(streamlines[2]) == false
    @test streamlines[1].escaped == true
    @test streamlines[2].escaped == false

    @testset "Mass loss rate" begin
        sl = streamlines[1]
        calculated = Qwind.compute_streamline_mdot(sl, bh.Rg)
        expected =
            4π *
            r[1] *
            widths[1] *
            bh.Rg^2 *
            n[1] *
            sqrt(vr1[1]^2 + vz1[1]^2) *
            C *
            M_P *
            0.61
        @test calculated ≈ expected rtol = 1e-2
        wind = Qwind.compute_wind_mdot(streamlines, bh.Rg)
        @test wind ≈ expected rtol = 1e-2 # 2nd streamline doesn't escape
    end

    @testset "Kinetic luminosity" begin
        sl = streamlines[1]
        calculated = Qwind.compute_kinetic_luminosity(sl, bh.Rg)
        vt = sqrt(0.3^2 + 0.2^2) * C
        expected = 0.5 * Qwind.compute_streamline_mdot(sl, bh.Rg) * vt^2
        @test calculated ≈ expected
        wind = Qwind.compute_kinetic_luminosity(streamlines, bh.Rg)
        @test wind ≈ expected # 2nd streamline doesn't escape
    end

    @testset "Wind momentum" begin
        sl = streamlines[1]
        calculated = Qwind.compute_momentum_rate(sl, bh.Rg)
        vt = sqrt(0.3^2 + 0.2^2) * C
        expected = Qwind.compute_streamline_mdot(sl, bh.Rg) * vt
        @test calculated ≈ expected
        wind = Qwind.compute_momentum_rate(streamlines, bh.Rg)
        @test wind ≈ expected # 2nd streamline doesn't escape
    end

    @testset "Max velocity" begin
        max_v = Qwind.compute_maximum_velocity(streamlines)
        @test max_v ≈ sqrt(0.3^2 + 0.2^2)
    end

    @testset "Wind properties" begin
        mu_n = 2
        props = Qwind.create_wind_properties(
            streamlines,
            bh,
            mu_nucleon = mu_n,
        )
        @test props["eddington_luminosity"] ≈
              compute_eddington_luminosity(bh)
        @test props["bolometric_luminosity"] ≈
              compute_bolometric_luminosity(bh)
        @test props["mass_accretion_rate"] ≈
              compute_mass_accretion_rate(bh)
        @test props["mass_loss"] ≈
              Qwind.compute_wind_mdot(streamlines, bh.Rg, mu_nucleon = mu_n)
        @test props["kinetic_luminosity"] ≈
              Qwind.compute_kinetic_luminosity(streamlines, bh.Rg, mu_nucleon = mu_n)
        @test props["max_velocity"] ≈ Qwind.compute_maximum_velocity(streamlines)
        @test props["momentum"] ≈
              Qwind.compute_momentum_rate(streamlines, bh.Rg, mu_nucleon = mu_n)
        @test props["mass_loss_fraction"] ≈
              props["mass_loss"] / props["mass_accretion_rate"]
    end

end


@testset "Test saving HDF5" begin

    model = Model(String(@__DIR__) * "/saver_tests_config.yaml")
    try
        rm("./saver_tests.hdf5", force = true)
    catch
    end
    iterations_dict = Dict(1 => Dict())
    integrators, streamlines =
        Qwind.run_integrators(model, iterations_dict, it_num = 1, parallel = false)
    trajectories = [Qwind.Trajectory(integ) for integ in integrators]
    Qwind.save_hdf5(integrators, streamlines, model, "./saver_tests.hdf5", 1)

    @testset "Save streamlines and trajectories" begin
        trajs_saved = load_trajectories("./saver_tests.hdf5")
        for (traj, traj_saved) in zip(trajectories, trajs_saved)
            for fieldname in fieldnames(Qwind.Trajectory)
                @test getfield(traj, fieldname) == getfield(traj_saved, fieldname)
            end
        end

        sl_saved = Streamlines("./saver_tests.hdf5")
        for (sl, sl_saved) in zip(streamlines, sl_saved)
            for fieldname in fieldnames(Qwind.Streamline)
                @test getfield(sl, fieldname) == getfield(sl_saved, fieldname)
            end
        end
    end

    @testset "Save wind hull" begin
        whull_saved = Hull("./saver_tests.hdf5")
        @test whull_saved.vertices == model.wind.wind_hull.vertices
    end

    @testset "Save density grid" begin
        recovered = DensityGrid("./saver_tests.hdf5")
        original = model.wind.density_grid
        for fieldname in fieldnames(DensityGrid)
            if fieldname == :iterator
                continue
            end
            @test getfield(original, fieldname) == getfield(recovered, fieldname)
        end
    end

    @testset "Save velocity grid" begin
        recovered = VelocityGrid("./saver_tests.hdf5")
        original = model.wind.velocity_grid
        for fieldname in fieldnames(VelocityGrid)
            if fieldname == :iterator
                continue
            end
            @test getfield(original, fieldname) == getfield(recovered, fieldname)
        end
    end

    #@testset "Save ionization grid" begin
    #    recovered = IonizationGrid("./saver_tests.hdf5")
    #    original = model.wind.ionization_grid
    #    for fieldname in fieldnames(IonizationGrid)
    #        if fieldname == :iterator
    #            continue
    #        end
    #        @test getfield(original, fieldname) == getfield(recovered, fieldname)
    #    end
    #end
end

