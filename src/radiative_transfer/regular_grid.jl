export RegularGrid, compute_uv_tau, compute_xray_tau

struct RegularGrid <: RadiativeTransfer
    radiation::Radiation
    kdtree::Union{KDTree,Nothing}
    density_interpolator::Union{DensityInterpolator,Nothing}
    vacuum_density::Float64
    atol::Float64
    rtol::Float64
    nr::Int
    nz::Int
    n_timesteps::Int
end

function RegularGrid(
    radiation::Radiation,
    integrators;
    vacuum_density = 1e2,
    atol = 1e-4,
    rtol = 1e-3,
    kernel_size = 1,
    nr = 500,
    nz = 501,
    n_timesteps = 10000,
)
    if integrators === nothing
        kdtree = nothing
    else
        kdtree = create_wind_kdtree(integrators, n_timesteps)
    end
    smoothed_grid = SmoothedGrid(kdtree, nr = nr, nz = nz, kernel_size = kernel_size)
    return RegularGrid(
        radiation,
        kdtree,
        smoothed_grid,
        vacuum_density,
        atol,
        rtol,
        nr,
        nz,
        n_timesteps,
    )
end

function RegularGrid(radiation::Radiation, config::Dict)
    rtc = config["radiative_transfer"]
    return RegularGrid(
        radiation,
        nothing,
        vacuum_density = rtc["vacuum_density"],
        atol = rtc["atol"],
        rtol = rtc["rtol"],
        kernel_size = rtc["kernel_size"],
        nr = rtc["nr"],
        nz = rtc["nz"],
        n_timesteps = rtc["n_timesteps"],
    )
end

get_density(regular_grid::RegularGrid, r, z) =
    get_density(regular_grid.density_interpolator, r, z)

function update_radiative_transfer(rt::RegularGrid, integrators)
    @info "Updating radiative transfer... "
    return RegularGrid(
        rt.radiation,
        integrators,
        vacuum_density = rt.vacuum_density,
        atol = rt.atol,
        rtol = rt.rtol,
        nr = rt.nr,
        nz = rt.nz,
        n_timesteps = rt.n_timesteps,
    )
end

compute_xray_tau(regular_grid::RegularGrid, ri, zi, rf, zf, xray_luminosity, Rg) =
    compute_xray_tau(regular_grid.density_interpolator, ri, zi, rf, zf, xray_luminosity, Rg)

compute_xray_tau(regular_grid::RegularGrid, r, z) = compute_xray_tau(
    regular_grid::RegularGrid,
    0.0,
    0.0,
    r,
    z,
    regular_grid.radiation.xray_luminosity,
    regular_grid.radiation.Rg,
)

compute_uv_tau(regular_grid::RegularGrid, ri, zi, rf, zf, Rg) =
    compute_uv_tau(regular_grid.density_interpolator, ri, zi, rf, zf, Rg)

compute_uv_tau(regular_grid::RegularGrid, rd, r, z) = compute_uv_tau(
    regular_grid::RegularGrid,
    rd, 
    0,
    r,
    z,
    regular_grid.radiation.Rg,
)
