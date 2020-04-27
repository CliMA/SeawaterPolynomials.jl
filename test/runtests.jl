using 
    Test,
    SeawaterPolynomials,
    SeawaterPolynomials.SecondOrderSeawaterPolynomials,
    SeawaterPolynomials.TEOS10

using SeawaterPolynomials: AbstractSeawaterPolynomial, BoussinesqEquationOfState

""" Test instantiation of a RoquetSeawaterPolynomial."""
function instantiate_roquet_polynomial(FT, coefficient_set)
    polynomial = RoquetSeawaterPolynomial(FT, coefficient_set)
    return typeof(polynomial) <: AbstractSeawaterPolynomial
end

""" Test instantiation of a RoquetEquationOfState."""
function instantiate_roquet_equation_of_state(FT, coefficient_set)
    eos = RoquetEquationOfState(FT, coefficient_set)
    return typeof(eos) <: BoussinesqEquationOfState
end

""" Test instantiation of TEOS10SeawaterPolynomial."""
function instantiate_teos10_polynomial(FT)
    polynomial = TEOS10SeawaterPolynomial(FT)
    return typeof(polynomial) <: AbstractSeawaterPolynomial
end

""" Test instantiation of a RoquetEquationOfState."""
function instantiate_teos10_equation_of_state(FT)
    eos = TEOS10EquationOfState(FT)
    return typeof(eos) <: BoussinesqEquationOfState
end

@testset "Second-order seawater polynomials" begin
    for coefficient_set in (
                            :Linear, 
                            :Cabbeling,
                            :CabbelingThermobaricity,
                            :Freezing,
                            :SecondOrder
                           )
    
        for FT in (Float64, Float32)
            @test instantiate_roquet_polynomial(FT, coefficient_set)
            @test instantiate_roquet_equation_of_state(FT, coefficient_set)

            eos = RoquetEquationOfState(FT)

            @test SeawaterPolynomials.ρ′(0, 0, 0, eos) == 0
            @test SeawaterPolynomials.haline_contraction(0, 0, 0, eos) == eos.seawater_polynomial.R₁₀₀
            @test SeawaterPolynomials.thermal_expansion(0, 0, 0, eos) == eos.seawater_polynomial.R₀₁₀
        end
    end

end

@testset "TEOS-10 seawater polynomials" begin
    for FT in (Float64, Float32)
        @test instantiate_teos10_polynomial(FT)
        @test instantiate_teos10_equation_of_state(FT)

        eos = TEOS10EquationOfState(FT)
        # Test/check values from Roquet et al. (2014).
        Θ = 10   # [C]
        S = 30   # [g/kg]
        Z = 1e3  # [m]

        τ = SeawaterPolynomials.TEOS10.τ(Θ)
        s = SeawaterPolynomials.TEOS10.s(S)
        ζ = SeawaterPolynomials.TEOS10.ζ(Z)

        @test SeawaterPolynomials.TEOS10.r₀(ζ) ≈ 4.59763035
        @test SeawaterPolynomials.TEOS10.r′(τ, s, ζ) ≈ 1022.85377

        @test SeawaterPolynomials.TEOS10.ρ′(Θ, S, Z, eos) ≈ 1027.45140
        @test SeawaterPolynomials.TEOS10.thermal_expansion(Θ, S, Z, eos) ≈ 0.179646281
        @test SeawaterPolynomials.TEOS10.haline_contraction(Θ, S, Z, eos) ≈ 0.765555368
    end
end
