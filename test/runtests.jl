using Test, IntegerMathUtils

@testset "iroot" begin
    for T in (Int32, Int64, BigInt)
        @test iroot(T(100), 2)   == T(10)
        @test iroot(T(101), 2)   == T(10)
        @test iroot(T(99),  2)   == T(9)
        @test iroot(T(1000), 3)  == T(10)
        @test iroot(T(1001), 3)  == T(10)
        @test iroot(T(999),  3)  == T(9)
        @test iroot(T(10000), 4) == T(10)
        @test iroot(T(10001), 4) == T(10)
        @test iroot(T(9999),  4) == T(9)
    end
    @test iroot(big(23)^50, 50)     == big(23)
    @test iroot(big(23)^50 + 1, 50) == big(23)
    @test iroot(big(23)^50 - 1, 50) == big(22)
end

@testset "ispower" begin
    for T in (Int32, Int64, BigInt)
        @test ispower(T(100))    == true
        @test ispower(T(1))      == true
        @test ispower(T(0))      == true
        @test ispower(T(12))     == false
        @test ispower(T(2^30))   == true
        @test ispower(T(5)^10)   == true
        @test ispower(T(2^30)+1) == false
        @test ispower(T(5)^10+1) == false
    end
    @test ispower(big(5)^40)     == true
    @test ispower(big(5)^40 + 1) == false
    @test ispower(6 * big(5)^40) == false
end

@testset "find_exponent" begin
    for T in (Int32, Int64, BigInt)
        @test find_exponent(T(100))    === 2
        @test find_exponent(T(1))      === 1
        @test find_exponent(T(0))      === 1
        @test find_exponent(T(12))     === 1
        @test find_exponent(T(2^30))   === 30
        @test find_exponent(T(5)^10)   === 10
        @test find_exponent(T(2^30)+1) === 1
        @test find_exponent(T(5)^10+1) === 1
    end
    @test find_exponent(big(5)^40)     === 40
    @test find_exponent(big(5)^40 + 1) === 1
    @test find_exponent(6*big(5)^40 )  === 1
end

@testset "is_probably_prime" begin
    for T in (Int32, Int64, BigInt)
        @test is_probably_prime(T(2)^7-1)  == true
        @test is_probably_prime(T(2)^13-1) == true
        @test is_probably_prime(T(2)^19-1) == true
        @test is_probably_prime(T(2)^27-1) == false
        @test is_probably_prime(T(2)^23-1) == false
        @test is_probably_prime(T(2)^30-1) == false
    end
    @test is_probably_prime(big(2)^127-1) == true
    @test is_probably_prime(big(2)^128-1) == false
end

@testset "kroneker" begin
    @test kronecker(4, 5) == kronecker(4, -5) == 1
    @test kronecker(1, 0) == kronecker(-1, 0) == 1
    @test kronecker(-4, -5) == -1
    @test kronecker(4, 6) == kronecker(-4, 0) == 0
end
