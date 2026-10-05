using
    Test,
    Random,
    SeawaterPolynomials,
    SeawaterPolynomials.SecondOrderSeawaterPolynomials,
    SeawaterPolynomials.TEOS10

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

@testset "TEOS-10 temperature conversions vs GibbsSeaWater" begin
    pinned_points = [
      # (Sᴬ,   T,      p,    label)
        (35.0, 10.0,    0.0, "mid Atlantic surface"),
        (34.5, 20.0,  100.0, "equatorial Pacific"),
        (34.9,  2.0, 4000.0, "NE Atlantic deep"),
        (33.0, -1.0, 1500.0, "Southern Ocean"),
        ( 8.0,  5.0,    0.0, "Baltic Sea"),
        (35.0, 25.0,    0.0, "Caribbean"),
    ]

    @testset "pinned points: $label" for (Sᴬ, T, p, label) in pinned_points
        Θ = GibbsSeaWater.gsw_ct_from_t(Sᴬ, T, p)
        @test Θ_from_θ(Sᴬ, T)    ≈ GibbsSeaWater.gsw_ct_from_pt(Sᴬ, T)    rtol=1e-12
        @test θ_from_T(Sᴬ, T, p) ≈ GibbsSeaWater.gsw_pt0_from_t(Sᴬ, T, p) rtol=1e-12
        @test Θ_from_T(Sᴬ, T, p) ≈ GibbsSeaWater.gsw_ct_from_t(Sᴬ, T, p)  rtol=1e-12
        @test θ_from_Θ(Sᴬ, Θ)    ≈ GibbsSeaWater.gsw_pt_from_ct(Sᴬ, Θ)    rtol=1e-12
    end

    @testset "random sweep" begin
        Random.seed!(20260429)
        for _ in 1:200
            Sᴬ = 30 + 10 * rand()
            θ  = -2 + 32 * rand()
            T  = -2 + 32 * rand()
            p  = 6000 * rand()
            Θ  = GibbsSeaWater.gsw_ct_from_pt(Sᴬ, θ)

            @test Θ_from_θ(Sᴬ, θ)    ≈ GibbsSeaWater.gsw_ct_from_pt(Sᴬ, θ)    rtol=1e-12
            @test θ_from_T(Sᴬ, T, p) ≈ GibbsSeaWater.gsw_pt0_from_t(Sᴬ, T, p) rtol=1e-12
            @test Θ_from_T(Sᴬ, T, p) ≈ GibbsSeaWater.gsw_ct_from_t(Sᴬ, T, p)  rtol=1e-12
            @test θ_from_Θ(Sᴬ, Θ)    ≈ GibbsSeaWater.gsw_pt_from_ct(Sᴬ, Θ)    rtol=1e-12
        end
    end
end

@testset "TEOS-10 freezing temperature" begin
    function gsw_freezing_derivatives(Sᴬ, p, saturation_fraction)
        salinity_derivative, pressure_derivative = Ref{Float64}(), Ref{Float64}()
        GibbsSeaWater.gsw_ct_freezing_first_derivatives_poly(Sᴬ, p, saturation_fraction, salinity_derivative, pressure_derivative)
        return salinity_derivative[], pressure_derivative[] * 1e4 # K/Pa → K/dbar
    end

    function central_difference(f, x, h)
        return (f(x + h) - f(x - h)) / 2h
    end

    pinned_points = [
      # (Sᴬ,   p,      label)
        (35.0,    0.0, "standard ocean surface"),
        (34.5,  500.0, "ice shelf base"),
        (34.7, 2000.0, "deep grounding line"),
        ( 8.0,    0.0, "Baltic Sea"),
        ( 0.0,    0.0, "fresh water"),
    ]

    @testset "pinned points: $label" for (Sᴬ, p, label) in pinned_points
        for saturation_fraction in (0, 0.5, 1)
            @test Θᶠ(Sᴬ, p, saturation_fraction) ≈ GibbsSeaWater.gsw_ct_freezing_poly(Sᴬ, p, saturation_fraction) rtol=1e-12
        end
        @test Θᶠ(Sᴬ, p) == Θᶠ(Sᴬ, p, 1)
    end

    @testset "random sweep" begin
        Random.seed!(20260929)
        for _ in 1:200
            Sᴬ = 1 + 41 * rand()
            p  = 3000 * rand()
            saturation_fraction = rand()

            @test Θᶠ(Sᴬ, p, saturation_fraction) ≈ GibbsSeaWater.gsw_ct_freezing_poly(Sᴬ, p, saturation_fraction) rtol=1e-12

            # GSW-C's salinity derivative has the wrong sign on one dissolved-air term, so it is only used air-free
            expected_salinity_derivative, expected_pressure_derivative = gsw_freezing_derivatives(Sᴬ, p, 0.0)
            @test ∂Θᶠ∂Sᴬ(Sᴬ, p, 0) ≈ expected_salinity_derivative rtol=1e-12
            @test ∂Θᶠ∂p(Sᴬ, p) ≈ expected_pressure_derivative rtol=1e-12

            expected_salinity_derivative = central_difference(S -> Θᶠ(S, p, saturation_fraction), Sᴬ, 1e-3)
            expected_pressure_derivative = central_difference(p′ -> Θᶠ(Sᴬ, p′, saturation_fraction), p, 1e-1)
            @test ∂Θᶠ∂Sᴬ(Sᴬ, p, saturation_fraction) ≈ expected_salinity_derivative rtol=1e-8
            @test ∂Θᶠ∂p(Sᴬ, p) ≈ expected_pressure_derivative rtol=1e-8
        end
    end

    @testset "Float32" begin
        Θ = Θᶠ(34.5f0, 500f0)
        @test Θ isa Float32
        @test Θ ≈ GibbsSeaWater.gsw_ct_freezing_poly(34.5, 500.0, 1.0) rtol=1e-6
        @test ∂Θᶠ∂Sᴬ(34.5f0, 500f0) isa Float32
        @test ∂Θᶠ∂p(34.5f0, 500f0) isa Float32
        @test Θᶠ(35, 0) isa Float64
    end
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
