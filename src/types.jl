# general pieces of a model

export DensityInterpolator,
    GridInterpolator,
    FMNoInterp,
    FMInterp,
    Relativistic,
    NoRelativistic,
    ConstantUVFraction,
    RadialUVFraction,
    TauUVCalculationFlag,
    TauUVCenter,
    TauUVDisk,
    NoTauUV,
    Boost,
    Thomson

abstract type WindInterpolator{T<:AbstractFloat} end
abstract type InitialConditions{T<:AbstractFloat} end
abstract type Grid{T<:AbstractFloat} end
abstract type InterpolationGrid{T<:AbstractFloat} end

# Radiation flags
abstract type Flag end

# Force multiplier interpolation
abstract type FMInterpolationFlag <: Flag end
struct FMNoInterp <: FMInterpolationFlag end
struct FMInterp <: FMInterpolationFlag end

# type of X-Ray opacity
abstract type XRayOpacityFlag <: Flag end
struct Thomson <: XRayOpacityFlag end
struct Boost <: XRayOpacityFlag end

# Relativistic correction to luminosity
abstract type RelativisticFlag <: Flag end
struct Relativistic <: RelativisticFlag end
struct NoRelativistic <: RelativisticFlag end

# Options for disc integral
abstract type UVFractionFlag <: Flag end
struct NoUVFraction <: UVFractionFlag end
struct UVFraction <: UVFractionFlag end

# Method to compute Tau UV
abstract type TauUVCalculationFlag end
struct TauUVCenter <: TauUVCalculationFlag end
struct TauUVDisk <: TauUVCalculationFlag end
struct NoTauUV <: TauUVCalculationFlag end

# Radiative transfer
abstract type CellIterator{T<:AbstractFloat} end
