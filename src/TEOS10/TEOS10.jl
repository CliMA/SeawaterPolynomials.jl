module TEOS10

export
    TEOS10SeawaterPolynomial,
    TEOS10EquationOfState,
    Θ_from_θ,
    Θ_from_T,
    θ_from_T,
    θ_from_Θ,
    freezing_conservative_temperature,
    freezing_conservative_temperature_salinity_derivative,
    freezing_conservative_temperature_pressure_derivative

using SeawaterPolynomials: AbstractSeawaterPolynomial, BoussinesqEquationOfState

import SeawaterPolynomials: ρ, ρ′, thermal_sensitivity, haline_sensitivity, with_float_type

include("density_computation.jl")
include("temperature_conversions.jl")
include("freezing_temperature.jl")

end # module
