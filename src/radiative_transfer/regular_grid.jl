using Qwind
export RegularGrid, compute_uv_tau, compute_xray_tau

struct RegularGrid{T} <: RadiativeTransfer{T}
    radiation::Radiation{T}
    density_interpolator::DensityInterpolator{T}
end

function RegularGrid(
    radiation::Radiation,
    integrators,
    density_interpolator_type::DataType;
    kwargs...
)
    interpolator = density_interpolator_type(integrators; kwargs...)
    return new(
        radiation,
        interpolator,
    )
end

function RegularGrid(radiation::Radiation, config::Dict)
    rtc = config[:radiative_transfer]
    intc = rtc[:density_interpolator]
    intc_type = pop!(intc, :type)
    interpolator = getfield(Qwind, Symbol(intc_type))(nothing; intc...)
    return RegularGrid(
        radiation,
        interpolator,
    )
end

get_density(regular_grid::RegularGrid, r, z) =
    get_density(regular_grid.density_interpolator, r, z)

function update_radiative_transfer(rt::RegularGrid, integrators)
    @info "Updating radiative transfer... "
    new_interp = update_density_interpolator(rt.density_interpolator, integrators)
    return RegularGrid(rt.radiation, new_interp)
end

compute_xray_tau(regular_grid::RegularGrid, ri, zi, rf, zf, xray_luminosity, Rg) =
    compute_xray_tau(regular_grid.density_interpolator.grid, ri, zi, rf, zf, xray_luminosity, Rg)

compute_xray_tau(regular_grid::RegularGrid, r, z) = compute_xray_tau(
    regular_grid::RegularGrid,
    0.0,
    6.0,
    r,
    z,
    regular_grid.radiation.xray_luminosity,
    regular_grid.radiation.Rg,
)

compute_uv_tau(regular_grid::RegularGrid, ri, zi, rf, zf, Rg) =
    compute_uv_tau(regular_grid.density_interpolator.grid, ri, zi, rf, zf, Rg)

compute_uv_tau(regular_grid::RegularGrid, rd, r, z) = compute_uv_tau(
    regular_grid::RegularGrid,
    rd, 
    0.0,
    r,
    z,
    regular_grid.radiation.Rg,
)
