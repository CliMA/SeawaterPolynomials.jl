module TEOS10

export
    TEOS10SeawaterPolynomial,
    TEOS10EquationOfState,
    Θ_from_θᴾ,
    Θ_from_T,
    θᴾ_from_T,
    θᴾ_from_Θ

using SeawaterPolynomials: AbstractSeawaterPolynomial, BoussinesqEquationOfState

import SeawaterPolynomials: ρ, ρ′, thermal_sensitivity, haline_sensitivity, with_float_type

include("density_computation.jl")
include("temperature_conversions.jl")

end # module
