module SeawaterPolynomials

export EquationOfState, SecondOrderSeawaterPolynomial, RoquetLinearSeawaterPolynomial,
       ρ, ρ′

abstract type AbstractSeawaterPolynomial end

struct EquationOfState{P, FT}
      reference_density :: FT
    seawater_polynomial :: P
end

@Base.kwdef struct SecondOrderSeawaterPolynomial{FT} <: AbstractSeawaterPolynomial
    R₀₁₀ :: FT
    R₀₀₁ :: FT
    R₀₂₀ :: FT
    R₀₁₁ :: FT
    R₂₀₀ :: FT
    R₁₀₁ :: FT
    R₁₁₀ :: FT
end

RoquetLinearSeawaterPolynomial(FT=Float64) =
    SecondOrderSeawaterPolynomial{FT}(
                                      R₀₁₀ = -1.775e-1,
                                      R₀₀₁ = 7.718e-1,
                                      R₀₂₀ = 0,
                                      R₀₁₁ = 0,
                                      R₂₀₀ = 0,
                                      R₁₀₁ = 0,
                                      R₁₁₀ = 0
                                     )

const EOS = EquationOfState

#ρ′(Θ, Sᴬ, Z, eos::EOS{<:SecondOrderSeawaterPolynomial}) = 

end # module
