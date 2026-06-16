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
                            :SecondOrder,
                            :SimplestRealistic
                           )

        for FT in (Float64, Float32)
            @test instantiate_roquet_polynomial(FT, coefficient_set)
            @test instantiate_roquet_equation_of_state(FT, coefficient_set)

            eos = RoquetEquationOfState(FT)

            @test SeawaterPolynomials.ρ′(0, 0, 0, eos) == 0
            @test SeawaterPolynomials.haline_sensitivity(0, 0, 0, eos) ==
                    eos.seawater_polynomial.R₁₀₀
            @test SeawaterPolynomials.thermal_sensitivity(0, 0, 0, eos) ==
                    eos.seawater_polynomial.R₀₁₀
        end
    end

end

@testset "TEOS-10 seawater polynomials" begin
    for FT in (Float64, Float32)
        @test instantiate_teos10_polynomial(FT)
        @test instantiate_teos10_equation_of_state(FT)

        eos = TEOS10EquationOfState(FT)
        # Test/check values from Roquet et al. (2014).
        Θ = 10.0    # [C]
        S = 30.0    # [g/kg]
        Z = - 1e3 # [m]

        τ = SeawaterPolynomials.TEOS10.τ(Θ)
        s = SeawaterPolynomials.TEOS10.s(S)
        ζ = SeawaterPolynomials.TEOS10.ζ(Z)

        @test SeawaterPolynomials.TEOS10.r₀(ζ) ≈ 4.59763035
        @test SeawaterPolynomials.TEOS10.r′(τ, s, ζ) ≈ 1022.85377

        @test SeawaterPolynomials.ρ(Θ, S, Z, eos) ≈ 1027.45140
        
        @test SeawaterPolynomials.TEOS10.thermal_sensitivity(Θ, S, Z, eos) ≈ 0.179646281
        @test SeawaterPolynomials.TEOS10.haline_sensitivity(Θ, S, Z, eos) ≈ 0.765555368
    end

end

@testset "Freezing temperature" begin
    import SeawaterPolynomials.TEOS10: conservative_temperature_freezing, in_situ_temperature_freezing

    # Reference values, with the geopotential height Z = -p [m ≈ dbar]. Each tuple is
    # (Sᴬ, Z, saturation_fraction, CT_freezing, t_freezing).
    #
    # The CT_freezing column is gsw.CT_freezing_poly(Sᴬ, p, sf) from the GSW-Python toolbox
    # (gsw 3.6.x); our `conservative_temperature_freezing` reproduces it to machine precision.
    #
    # The t_freezing column is the standalone in-situ-temperature freezing polynomial (the
    # direct polynomial documented in GSW's gsw_t_freezing_poly source). It is accurate to
    # ~3e-4 °C against the exact gsw.t_freezing; the GSW gsw_t_freezing_poly function instead
    # routes through a CT→t conversion (requiring the full Gibbs function), so it differs from
    # this direct polynomial at that ~1e-4 °C level.
    reference_values = (
        (35.0,  -300.0, 0.0, -2.13822778628046,   -2.13578432707088),
        (35.0,  -300.0, 1.0, -2.14012815279748,   -2.13768667371971),
        (30.0,     0.0, 0.0, -1.61856867110253,   -1.62756358904631),
        (33.0,  -100.0, 0.5, -1.86841270639344,   -1.87230317368548),
        (34.7, -1000.0, 0.0, -2.67353983508005,   -2.65703117757126),
        ( 0.0,     0.0, 0.0,  0.0179470643279687,  0.002519),
    )

    for (Sᴬ, Z, sf, Θf, Tf) in reference_values
        @test conservative_temperature_freezing(Sᴬ, Z, saturation_fraction=sf) ≈ Θf
        @test in_situ_temperature_freezing(Sᴬ, Z, saturation_fraction=sf) ≈ Tf
    end

    # Default saturation_fraction is 0
    @test conservative_temperature_freezing(35.0, -300.0) == conservative_temperature_freezing(35.0, -300.0, saturation_fraction=0)
    @test in_situ_temperature_freezing(35.0, -300.0) == in_situ_temperature_freezing(35.0, -300.0, saturation_fraction=0)

    # Float type is preserved
    for FT in (Float64, Float32)
        @test conservative_temperature_freezing(FT(35), FT(-300)) isa FT
        @test in_situ_temperature_freezing(FT(35), FT(-300)) isa FT
    end
end

@testset "show" begin
    show_polynomial_string = repr(RoquetSeawaterPolynomial(:SecondOrder))
    R₀₁₀ = 0.182e-1
    R₁₀₀ = 8.078e-1
    R₀₂₀ = 4.937e-3
    R₀₁₁ = 2.4677e-5
    R₂₀₀ = 1.115e-4
    R₁₀₁ = 8.241e-6
    R₁₁₀ = 2.446e-3
    test_polynomial_string =
        "ρ' = $(eval(R₁₀₀)) Sᴬ + $(eval(R₀₁₀)) Θ - $(eval(R₀₂₀)) Θ² - $(eval(R₀₁₁)) Θ Z - $(eval(R₂₀₀)) Sᴬ² - $(eval(R₁₀₁)) Sᴬ Z - $(eval(R₁₁₀)) Sᴬ Θ"
    @test show_polynomial_string == test_polynomial_string
end

@testset "with_float_type" begin
    for (FT, FT2) in zip((Float32, Float64), (Float64, Float32))
        eos = TEOS10EquationOfState(FT)
        @test eltype(eos) == FT

        @show FT2
        eos = SeawaterPolynomials.with_float_type(FT2, eos)
        @test eltype(eos) == FT2

        for coefficient_set in (:Linear,
                                :Cabbeling,
                                :CabbelingThermobaricity,
                                :Freezing,
                                :SecondOrder,
                                :SimplestRealistic)

            eos = RoquetEquationOfState(FT, coefficient_set)
            @test eltype(eos) == FT

            eos = SeawaterPolynomials.with_float_type(FT2, eos)
            @test eltype(eos) == FT2
        end
    end
end
