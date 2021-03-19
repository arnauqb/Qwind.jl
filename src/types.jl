# general pieces of a model

export DensityInterpolator, GridInterpolator, FMNoInterp, FMInterp, Relativistic, NoRelativistic

abstract type Radiation{T<:AbstractFloat} end
abstract type RadiativeTransfer{T<:AbstractFloat} end
abstract type Interpolator{T<:AbstractFloat} end
abstract type InitialConditions{T<:AbstractFloat} end
abstract type Grid{T<:AbstractFloat} end
abstract type InterpolationGrid{T<:AbstractFloat} end

# Radiation
# Force multiplier interpolation
abstract type Flag end
abstract type FMInterpolationType <: Flag end
struct FMNoInterp <: FMInterpolationType end
struct FMInterp <: FMInterpolationType end
# Relativistic correction to luminosity
abstract type FluxCorrection <: Flag end
struct Relativistic <: FluxCorrection end
struct NoRelativistic <: FluxCorrection end

# Radiative transfer
abstract type CellIterator{T<:AbstractFloat} end
#abstract type LogGridIterator <: CellIterator end
