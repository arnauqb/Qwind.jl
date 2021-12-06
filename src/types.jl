# general pieces of a model

export FMInterpolationFlag,
    FMNoInterp,
    FMInterp,
    RelativisticFlag,
    Relativistic,
    NoRelativistic,
    ConstantUVFraction,
    RadialUVFraction,
    TauUVCalculationFlag,
    TauUVCenter,
    TauUVDisk,
    NoTauUV,
    OpacityFlag,
    BoostOpacity,
    ThomsonOpacity,
    NoOpacity,
    XRayScatteringFlag,
    Scattering,
    NoScattering,
    ICFlag,
    CAKMode,
    UniformMode,
    UpdateGridFlag,
    AverageGrid,
    ReplaceGrid


abstract type InitialConditions{T<:AbstractFloat} end
abstract type InterpolationGrid end

# Radiation flags
abstract type Flag end

# Force multiplier interpolation
abstract type FMInterpolationFlag <: Flag end
struct FMNoInterp <: FMInterpolationFlag end
struct FMInterp <: FMInterpolationFlag end

# opacity types
abstract type OpacityFlag <: Flag end
struct NoOpacity <: OpacityFlag end
struct ThomsonOpacity <: OpacityFlag end
struct BoostOpacity <: OpacityFlag end
struct BoostOpacity2 <: OpacityFlag end

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

# Method to update grids
abstract type UpdateGridFlag end
struct AverageGrid <: UpdateGridFlag end
struct ReplaceGrid <: UpdateGridFlag end

# X-ray scattering
abstract type XRayScatteringFlag end
struct Scattering <: XRayScatteringFlag end
struct NoScattering <: XRayScatteringFlag end

# Initial conditions flags
abstract type ICFlag end
struct CAKMode <: ICFlag end
struct UniformMode <: ICFlag end
