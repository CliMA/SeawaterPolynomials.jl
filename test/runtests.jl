using 
    Test,
    SeawaterPolynomials,
    SeawaterPolynomials.SecondOrderSeawaterPolynomials

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
