module TEOS10

export
    TEOS10SeawaterPolynomial,
    TEOS10EquationOfState,
    Θ_from_θᴾ,
    Θ_from_T,
    θᴾ_from_T,
    Sᴬ_from_Sᴾ

using LazyArtifacts

using SeawaterPolynomials: AbstractSeawaterPolynomial, BoussinesqEquationOfState

import SeawaterPolynomials: ρ, ρ′, thermal_sensitivity, haline_sensitivity, with_float_type

include("density_computation.jl")
include("temperature_conversions.jl")
include("salinity_conversions.jl")

# Populate the SAAR reference atlas from the bundled binary blob at module load time. The Ref holds an
# immutable `SAARAtlas`, set once here and read-only thereafter.
function __init__()
    SAAR_ATLAS[] = load_saar_atlas(saar_data_path())
    return nothing
end

end # module
