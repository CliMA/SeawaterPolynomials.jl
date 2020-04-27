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
        end
    end
end

@testset "TEOS-10 seawater polynomials" begin
    for FT in (Float64, Float32)
        @test instantiate_teos10_polynomial(FT)
        @test instantiate_teos10_equation_of_state(FT)
    end
end
