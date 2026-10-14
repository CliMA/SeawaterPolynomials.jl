module TEOS10

export
    TEOS10SeawaterPolynomial,
    TEOS10EquationOfState,
    Θ_from_θ,
    Θ_from_T,
    θ_from_T,
    θ_from_Θ,
    Sᴬ_from_Sᴾ

using Artifacts

import Adapt

using SeawaterPolynomials: AbstractSeawaterPolynomial, BoussinesqEquationOfState

import SeawaterPolynomials: ρ, ρ′, thermal_sensitivity, haline_sensitivity, with_float_type

include("density_computation.jl")
include("temperature_conversions.jl")
include("salinity_conversions.jl")

function __init__()
    SAAR_ATLAS[] = SAARAtlas(joinpath(artifact"gsw_saar_data", "gsw_saar_data.bin"))
    return nothing
end

end # module
