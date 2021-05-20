using Qwind
export RegularGrid, compute_uv_tau, compute_xray_tau

struct RegularGrid{T} <: RadiativeTransfer{T}
    radiation::Radiation{T}
    interpolator::Interpolator{T}
end

function RegularGrid(
    radiation::Radiation,
    integrators,
    interpolator_type::DataType;
    kwargs...
)
    interpolator = interpolator_type(integrators; kwargs...)
    return new(
        radiation,
        interpolator,
    )
end

function RegularGrid(radiation::Radiation, config::Dict)
    rtc = config[:radiative_transfer]
    intc = rtc[:interpolator]
    intc_type = pop!(intc, :type)
    interpolator = getfield(Qwind, Symbol(intc_type))(nothing; intc...)
    return RegularGrid(
        radiation,
        interpolator,
    )
end

get_density(regular_grid::RegularGrid, r, z) =
    get_density(regular_grid.interpolator, r, z)

function update_radiative_transfer(rt::RegularGrid, integrators)
    @info "Updating radiative transfer... "
    flush()
    new_interp = update_interpolator(rt.interpolator, integrators)
    return RegularGrid(rt.radiation, new_interp)
end

compute_xray_tau(regular_grid::RegularGrid, ri, zi, rf, zf, xray_luminosity, Rg) =
    compute_xray_tau(regular_grid.interpolator.density_grid, ri, zi, rf, zf, xray_luminosity, Rg)

compute_xray_tau(regular_grid::RegularGrid, z_xray, r, z) = compute_xray_tau(
    regular_grid::RegularGrid,
    0.0,
    z_xray,
    r,
    z,
    regular_grid.radiation.xray_luminosity,
    regular_grid.radiation.Rg,
)

compute_uv_tau(regular_grid::RegularGrid, ri, zi, rf, zf, Rg) =
    compute_uv_tau(regular_grid.interpolator.density_grid, ri, zi, rf, zf, Rg)

compute_uv_tau(regular_grid::RegularGrid, rd, r, z) = compute_uv_tau(
    regular_grid::RegularGrid,
    rd, 
    regular_grid.radiation.zh,
    r,
    z,
    regular_grid.radiation.Rg,
)
