using Test
using Random
using SeawaterPolynomials
using SeawaterPolynomials.SecondOrderSeawaterPolynomials
using SeawaterPolynomials.TEOS10

using SeawaterPolynomials: AbstractSeawaterPolynomial, BoussinesqEquationOfState

import GibbsSeaWater

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

@testset "TEOS-10 conversions vs GibbsSeaWater" begin
    pinned_points = [
      # (Sᴬ,    T,       p,   Sᴾ,     λ,     φ,    label)
        (35.0,  10.0,    0.0, 35.0,  -30.0,  45.0, "mid Atlantic surface"),
        (34.5,  20.0,  100.0, 34.5,  200.0,  -5.0, "equatorial Pacific"),
        (34.9,   2.0, 4000.0, 34.9,   20.0,  60.0, "NE Atlantic deep"),
        (33.0,  -1.0, 1500.0, 33.0,  -20.0, -60.0, "Southern Ocean"),
        ( 8.0,   5.0,    0.0,  8.0,   20.0,  60.0, "Baltic Sea"),
        (35.0,  25.0,    0.0, 35.0,  -78.0,   8.0, "Caribbean / Panama region"),
        (34.8,  15.0,  100.0, 34.8,  -82.0,  10.0, "north of Panama isthmus"),
    ]

    @testset "pinned points: $label" for (Sᴬ, T, p, Sᴾ, λ, φ, label) in pinned_points
        @test Θ_from_θᴾ(Sᴬ, T)        ≈ GibbsSeaWater.gsw_ct_from_pt(Sᴬ, T)        atol=0 rtol=1e-12
        @test θᴾ_from_T(Sᴬ, T, p)     ≈ GibbsSeaWater.gsw_pt0_from_t(Sᴬ, T, p)     atol=0 rtol=1e-12
        @test Θ_from_T(Sᴬ, T, p)      ≈ GibbsSeaWater.gsw_ct_from_t(Sᴬ, T, p)      atol=0 rtol=1e-12
        @test Sᴬ_from_Sᴾ(Sᴾ, p, λ, φ) ≈ GibbsSeaWater.gsw_sa_from_sp(Sᴾ, p, λ, φ)  atol=0 rtol=1e-12
    end

    # Random sweep over the standard oceanographic envelope. 
    # The fixed seed keeps CI deterministic.
    @testset "random sweep" begin
        Random.seed!(20260429)
        for _ in 1:200
            Sᴬ = 30 + 10 * rand()
            θᴾ = -2 + 32 * rand()
            T  = -2 + 32 * rand()
            p  = 6000 * rand()
            Sᴾ = 30 + 10 * rand()
            λ  = 360 * rand()
            φ  = -85 + 175 * rand()

            @test Θ_from_θᴾ(Sᴬ, θᴾ)   ≈ GibbsSeaWater.gsw_ct_from_pt(Sᴬ, θᴾ)   atol=0 rtol=1e-12
            @test θᴾ_from_T(Sᴬ, T, p) ≈ GibbsSeaWater.gsw_pt0_from_t(Sᴬ, T, p) atol=0 rtol=1e-12
            @test Θ_from_T(Sᴬ, T, p)  ≈ GibbsSeaWater.gsw_ct_from_t(Sᴬ, T, p)  atol=0 rtol=1e-12

            Sᴬˢ = Sᴬ_from_Sᴾ(Sᴾ, p, λ, φ)
            Sᴬᴳ = GibbsSeaWater.gsw_sa_from_sp(Sᴾ, p, λ, φ)

            if !isnan(Sᴬˢ) && !isnan(Sᴬᴳ)
                @test Sᴬˢ ≈ Sᴬᴳ atol=0 rtol=1e-12
            end
        end
    end
end

@testset "convert" begin
    for (FT, FT2) in zip((Float32, Float64), (Float64, Float32))
        teos10 = instantiate_teos10_polynomial(FT)
        test10 = convert(FT2, teos10)
        @test eltype(teos10) == FT2

        for coefficient_set in (:Linear,
                                :Cabbeling,
                                :CabbelingThermobaricity,
                                :Freezing,
                                :SecondOrder,
                                :SimplestRealistic)

            eos = instantiate_roquet_equation_of_state(FT, coefficient_set)
            eos = convert(FT2, teos10)
            @test eltype(teos10) == FT2
        end
    end
end