using HDF5, Printf
import ConcaveHull
export IonizationGrid, get_ionization_parameter, interpolate_ionization_parameter

#struct IonizationGrid{T,U,V} <: InterpolationGrid
#    r_range::Vector{T}
#    z_range::Vector{T}
#    grid::Array{T,2}
#    nr::Union{U,String}
#    nz::U
#    iterator::GridIterator{T,U,V}
#    interpolator::Any
#    function IonizationGrid(r_range, z_range, grid, nr = nothing, nz = nothing)
#        if nr === nothing
#            nr = length(r_range)
#        end
#        if nz === nothing
#            nz = length(z_range)
#        end
#        iterator = GridIterator(r_range, z_range)
#        interpolator =
#            Interpolations.interpolate((r_range, z_range), grid, Gridded(Linear()))
#        interpolator = Interpolations.extrapolate(interpolator, 1e15)
#        return new{typeof(r_range[1]),Int,Bool}(
#            r_range,
#            z_range,
#            grid,
#            nr,
#            nz,
#            iterator,
#            interpolator,
#        )
#    end
#end
#
#IonizationGrid(grid_data::Dict) = IonizationGrid(
#    grid_data["r"],
#    grid_data["z"],
#    grid_data["grid"],
#    grid_data["nr"],
#    grid_data["nz"],
#)
#
#function IonizationGrid(h5_path::String, it_num)
#    it_name = @sprintf "iteration_%03d" it_num
#    grid_data = h5open(h5_path, "r") do file
#        read(file, it_name * "/ionization_grid")
#    end
#    return IonizationGrid(grid_data)
#end
#
#function IonizationGrid(h5_path::String)
#    it_keys = h5open(h5_path, "r") do file
#        keys(read(file))
#    end
#    it_nums = [parse(Int, split(key, "_")[end]) for key in it_keys]
#    return IonizationGrid(h5_path, maximum(it_nums))
#end
#
#function IonizationGrid(nr::Union{String,Int}, nz::Int, vaccum_ion::Float64)
#    r_range = [-1.0, 0.0]
#    z_range = [-1.0, 0.0]
#    grid = vaccum_ion .* [[1.0, 1.0] [1.0, 1.0]]
#    return IonizationGrid(r_range, z_range, grid, nr, nz)
#end
#
#function compute_change_in_grid(old_grid, new_grid)
#    return mean(abs.(old_grid .- new_grid) ./ old_grid)
#end
#
#function compute_initial_ionization_grid(density_grid, source_luminosity, Rg, z_xray)
#    rr = density_grid.r_range[1:end-1]
#    zz = density_grid.z_range[1:end-1]
#    ret = zeros(length(rr), length(zz))
#    for (i, r) in enumerate(rr)
#        for (j, z) in enumerate(zz)
#            dist2 = (r^2 + (z - z_xray)^2) * Rg^2
#            tau = compute_optical_depth(
#                density_grid.iterator,
#                density_grid,
#                IonizationGrid(density_grid.nr, density_grid.nz, 1e10),
#                ThomsonOpacity(),
#                ri = 0,
#                phii = 0,
#                zi = z_xray,
#                rf = r,
#                phif = 0.0,
#                zf = z,
#                Rg=Rg
#            )
#            density = density_grid.grid[i, j] 
#            ret[i, j] = source_luminosity * exp(-tau) / density / dist2
#        end
#    end
#    ret = ret ./ density_grid.grid[1:(end - 1), 1:(end - 1)]
#    return IonizationGrid(rr, zz, ret)
#end
#
#function IonizationGrid(
#    density_grid::DensityGrid;
#    xray_luminosity,
#    Rg,
#    z_xray,
#    mu_nucleon = 0.61,
#    mu_electron = 1.13,
#    absorption_opacity=BoostOpacity(),
#    include_scattered=true,
#    parallel=true,
#)
#    @info "Updating ionization grid..."
#    xi0 = compute_initial_ionization_grid(density_grid, xray_luminosity, Rg, z_xray)
#    xi_grid = zeros(xi0.nr, xi0.nz)
#    for k in 1:10
#        @info "Iteration $k..."
#        total_flux = compute_total_flux_grid(
#            density_grid,
#            xi0,
#            parallel = parallel,
#            Rg = Rg,
#            mu_nucleon = mu_nucleon,
#            mu_electron = mu_electron,
#            z_xray = z_xray,
#            xray_luminosity = xray_luminosity,
#            absorption_opacity=absorption_opacity,
#            include_scattered=include_scattered
#        )
#        change = compute_change_in_grid(xi0.grid, xi_grid)
#        xi_grid = total_flux ./ density_grid.grid[1:(end - 1), 1:(end - 1)]
#        if abs(1 - change) < 1e-1
#            break
#        end
#    end
#    return IonizationGrid(
#        density_grid.r_range[1:(end - 1)],
#        density_grid.z_range[1:(end - 1)],
#        xi_grid,
#    )
#end
#
#
#function get_ionization_parameter(grid::IonizationGrid, r, z)
#    ridx = searchsorted_nearest(grid.r_range, r)
#    zidx = searchsorted_nearest(grid.z_range, z)
#    return grid.grid[ridx, zidx]
#end
#
#get_ionization_parameter(grid::IonizationGrid, point::Vector{Float64}) =
#    get_ionization_parameter(grid, point[1], point[2])
#
#function interpolate_ionization_parameter(grid::IonizationGrid, r, z)
#    return grid.interpolator(r, z)
#end
