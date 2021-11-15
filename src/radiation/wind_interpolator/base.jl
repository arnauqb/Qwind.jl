using PyCall
import ConcaveHull, Interpolations, Sundials
export WindInterpolator, get_density

struct WindInterpolator{T,U,V}
    wind_hull::ConcaveHull.Hull
    density_grid::DensityGrid{T,U,V}
    velocity_grid::VelocityGrid{T,U,V}
    #ionization_grid::IonizationGrid{T,U,V}
    vacuum_density::T
    update_grid_flag::UpdateGridFlag
end

WindInterpolator(nr, nz; vacuum_density = 1e2, update_grid_flag = AverageGrid()) =
    WindInterpolator(
        Hull(),
        DensityGrid(nr, nz, vacuum_density),
        VelocityGrid(nr, nz, 0.0),
        #IonizationGrid(nr, nz, 1e15),
        vacuum_density,
        update_grid_flag,
    )

function WindInterpolator(config::Dict)
    update_method = get(config, :update_grid_method, "average")
    if update_method == "average"
        update_method = AverageGrid()
    elseif update_method == "replace"
        update_method = ReplaceGrid()
    else
        error("Grid update method $update_method not recoginsed")
    end
    return WindInterpolator(
        config[:nr],
        config[:nz],
        vacuum_density = config[:vacuum_density],
        update_grid_flag = update_method,
    )
end

function WindInterpolator(
    streamlines::Union{Streamlines,Nothing};
    nr = "auto",
    nz = 50,
    vacuum_density = 1e2,
    xray_luminosity,
    Rg,
    include_scattered_flux=true,
    mu_nucleon = 0.61,
    mu_electron = 1.17,
    z_xray,
    update_grid_flag = Average(),
)
    if nr != "auto"
        nr = Int(nr)
    end
    nz = Int(nz)
    if streamlines === nothing
        hull = Hull()
        density_grid = DensityGrid(nr, nz, vacuum_density)
        velocity_grid = VelocityGrid(nr, nz, 0.0)
        #ionization_grid = IonizationGrid(nr, nz, 1e15)
    else
        hull = Hull(streamlines)
        density_grid = DensityGrid(streamlines, hull, nr = nr, nz = nz)
        velocity_grid = VelocityGrid(streamlines, hull, nr = nr, nz = nz)
        #ionization_grid = IonizationGrid(
        #    density_grid,
        #    xray_luminosity = xray_luminosity,
        #    Rg = Rg,
        #    mu_nucleon = mu_nucleon,
        #    mu_electron = mu_electron,
        #    z_xray = z_xray,
        #    include_scattered = include_scattered_flux,
        #)
    end
    return WindInterpolator(
        hull,
        density_grid,
        velocity_grid,
        #ionization_grid,
        vacuum_density,
        update_grid_flag,
    )
end


function update_wind_interpolator(wi::WindInterpolator, dgrid::DensityGrid)
    return WindInterpolator(
        wi.wind_hull,
        dgrid,
        wi.velocity_grid,
        #wi.ionization_grid,
        wi.vacuum_density,
        wi.update_grid_flag,
    )
end

function update_wind_interpolator(wi::WindInterpolator, vgrid::VelocityGrid)
    return WindInterpolator(
        wi.wind_hull,
        wi.density_grid,
        vgrid,
        #wi.ionization_grid,
        wi.vacuum_density,
        wi.update_grid_flag,
    )
end

#function update_wind_interpolator(wi::WindInterpolator, xigrid::IonizationGrid)
#    return WindInterpolator(
#        wi.wind_hull,
#        wi.density_grid,
#        wi.velocity_grid,
#        xigrid,
#        wi.vacuum_density,
#        wi.update_grid_flag,
#    )
#end

function get_density(wi::WindInterpolator, r, z)
    if is_point_in_wind(wi, [r, z])
        return wi.grid.interpolator(r, z)
    else
        return wi.vacuum_density
    end
end
get_density(wi::WindInterpolator, point) = get_density(wi, point[1], point[2])

function GridIterator(interpolator::WindInterpolator, ri, zi, rf, zf)
    return GridIterator(
        interpolator.density_grid.r_range,
        interpolator.density_grid.z_range,
        ri,
        zi,
        rf,
        zf,
    )
end
