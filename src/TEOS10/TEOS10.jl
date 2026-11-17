module TEOS10

export
    TEOS10SeawaterPolynomial,
    TEOS10EquationOfState,
    θ_from_θᴾ,
    θ_from_T,
    θᴾ_from_T,
    θᴾ_from_θ,
    θᶠʳ,
    ∂θᶠʳ∂Sᴬ,
    ∂θᶠʳ∂p

using SeawaterPolynomials: AbstractSeawaterPolynomial, BoussinesqEquationOfState

import SeawaterPolynomials: ρ, ρ′, thermal_sensitivity, haline_sensitivity, with_float_type

include("density_computation.jl")
include("temperature_conversions.jl")
include("freezing_temperature.jl")

end # module
