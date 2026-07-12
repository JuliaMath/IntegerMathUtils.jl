using Test, IntegerMathUtils
using Random

@testset "IntegerMathUtils" begin

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
        @test iroot(T(-8), 3) == T(-2)
        @test_throws DomainError iroot(T(-8), 4)
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

@testset "iroot(x, 2) matches Base.isqrt" begin
    for T in (Int32, UInt32, Int64, UInt64, Int128, UInt128)
        @test iroot(typemax(T), 2) == Base.isqrt(typemax(T))
    end
    # perfect squares and neighbors, including near the Float64 exactness edge (2^52)
    for a in (94906263, 94906264, 94906265, 94906266, 3037000499, 2^26, 2^26 - 1)
        @test iroot(Int64(a)^2, 2)     == a
        @test iroot(Int64(a)^2 + 1, 2) == a
        @test iroot(Int64(a)^2 - 1, 2) == a - 1
    end
    rng = MersenneTwister(42)
    for T in (UInt32, Int64, UInt64, Int128, UInt128), _ in 1:200
        x = rand(rng, T) & typemax(T)  # nonnegative
        @test iroot(x, 2) === T(Base.isqrt(x))
    end
end

@testset "iroot extended semantics" begin
    for T in (Int32, Int64, Int128, BigInt)
        @test iroot(T(0), 2) == 0
        @test iroot(T(0), 3) == 0
        @test iroot(T(1), 100) == 1
        @test iroot(T(-9), 3) == -2   # truncated toward zero, like GMP
        @test iroot(T(-1), 5) == -1
        @test_throws DomainError iroot(T(5), 0)
        @test_throws DomainError iroot(T(5), -2)
        @test_throws DomainError iroot(T(-4), 2)
    end
    @test iroot(typemin(Int64), 63) == -2
    @test iroot(typemin(Int8), 7) == -2
    @test iroot(typemax(UInt128), 2) == UInt128(2)^64 - 1
    @test iroot(Int128(2)^126, 63) == 4
    @test iroot(typemax(UInt128), 127) == 2
    @test iroot(Int8(127), 100) == 1
end

@testset "rootrem" begin
    for T in (Int32, Int64, Int128, BigInt)
        r = rootrem(T(1001), 3)
        @test r == (10, 1)
        @test r isa Tuple{T,T}
        @test rootrem(T(1000), 3) == (10, 0)
        @test rootrem(T(-1001), 3) == (-10, -1)
        @test rootrem(T(-1000), 3) == (-10, 0)
        @test_throws DomainError rootrem(T(-4), 2)
        @test_throws DomainError rootrem(T(5), 0)
    end
    @test rootrem(typemin(Int64), 63) == (-2, 0)
end

@testset "negative perfect powers" begin
    for T in (Int32, Int64, Int128, BigInt)
        @test ispower(T(-8))  == true
        @test ispower(T(-4))  == false
        @test ispower(T(-1))  == true
        @test ispower(T(-16)) == false   # even exponents don't count for negatives
        @test ispower(T(-64)) == true    # (-4)^3
        @test find_exponent(T(-8))  === 3
        @test find_exponent(T(-4))  === 1
        @test find_exponent(T(-1))  === 1
        @test find_exponent(T(-64)) === 3
    end
    @test ispower(typemin(Int64)) == true       # (-2)^63
    @test find_exponent(typemin(Int64)) === 63
    @test find_exponent(typemin(Int8))  === 7
    @test find_exponent(-big(3)^15) === 15
end

@testset "perfect powers near typemax" begin
    @test ispower(Int64(3)^39) == true
    @test find_exponent(Int64(3)^39) === 39
    @test ispower(UInt64(3)^40) == true
    @test find_exponent(UInt64(3)^40) === 40
    @test ispower(Int64(2)^62) == true
    @test find_exponent(Int64(2)^62) === 62
    for x in (Int64(3)^39 + 1, Int64(3)^39 - 1, Int64(2)^62 + 1)
        @test ispower(x) == ispower(big(x))
        @test find_exponent(x) == find_exponent(big(x))
    end
end

@testset "cross-check native against BigInt/GMP" begin
    rng = MersenneTwister(1234)
    types = (Int8, UInt8, Int16, UInt16, Int32, UInt32, Int64, UInt64, Int128, UInt128)
    for T in types, _ in 1:200
        x = rand(rng, T)
        for n in (1, 2, 3, 5, 7, rand(rng, 2:130))
            if x < 0 && iseven(n)
                @test_throws DomainError iroot(x, n)
                @test_throws DomainError rootrem(x, n)
            else
                a = iroot(x, n)
                @test a isa T
                mx, ma = abs(big(x)), abs(big(a))
                @test ma^n <= mx < (ma + 1)^n
                @test iroot(big(x), n) == big(a)
                root, rem = rootrem(x, n)
                @test root === a
                @test rem isa T
                @test big(rem) == big(x) - big(a)^n
            end
        end
        @test ispower(x) == ispower(big(x))
        e = find_exponent(x)
        if abs(big(x)) > 1
            @test (e > 1) == ispower(big(x))
            @test big(iroot(x, e))^e == big(x)
        end
    end
end

@testset "overflow-check band (root == 2 near typemax)" begin
    # These are the cases that prove _pow_leq's overflow check is needed:
    # 2^n <= typemax < 3^n, so the true root is 2 while (root+1)^n overflows.
    @test iroot(typemax(UInt8), 7) == 2       # 3^7 == 139 (mod 256) <= 255
    @test iroot(Int8(127), 6) == 2
    for n in 41:63
        @test iroot(typemax(UInt64), n) == UInt64(2)
        @test iroot(typemax(UInt64), n) == UInt64(iroot(big(typemax(UInt64)), n))
    end
    for n in 85:127
        @test iroot(typemax(UInt128), n) == UInt128(2)
    end
end

@testset "near-integer reject: bᵏ and neighbors near typemax" begin
    for T in (Int32, UInt32, Int64, UInt64, Int128, UInt128)
        hi = big(typemax(T))
        for k in (3, 5, 7, 11, 3, 13)
            b = Int(floor(hi^(1/k)))
            for bb in (b-1, b, b+1)
                bb < 2 && continue
                p = big(bb)^k
                p > hi && continue
                x = T(p)
                @test ispower(x) == ispower(p)
                @test find_exponent(x) == find_exponent(p)
                # exact power must be detected, and near-misses rejected
                @test ispower(x) == true
                if x < typemax(T)
                    @test ispower(x + one(T)) == ispower(p + 1)
                end
                if x > T(2)
                    @test ispower(x - one(T)) == ispower(p - 1)
                end
            end
        end
    end
end

@testset "generic (non-BitInteger) native root path" begin
    # Public iroot(::BigInt) uses GMP, but the internal _iroot_nonneg / _iroot_approx
    # run the generic T<:Integer code (Newton + rescaled seed / UInt128 shortcut).
    # Exercise them on BigInt and cross-check against GMP.
    IMU = IntegerMathUtils
    rng = MersenneTwister(2024)
    for _ in 1:300
        x = rand(rng, big(0):(big(2)^rand(rng, 1:400)))
        for n in (2, 3, 5, 7, rand(rng, 2:80))
            @test IMU._iroot_nonneg(x, n) == iroot(x, n)
        end
    end
    # explicitly hit both branches of the generic _iroot_approx and their edges
    @test IMU._iroot_nonneg(big(2)^127, 3)   == iroot(big(2)^127, 3)    # <= typemax(UInt128)
    @test IMU._iroot_nonneg(big(2)^128, 3)   == iroot(big(2)^128, 3)    # just over
    @test IMU._iroot_nonneg(big(10)^100, 5)  == iroot(big(10)^100, 5)   # large, rescaled seed
    @test IMU._iroot_nonneg(big(2)^1000, 3)  == iroot(big(2)^1000, 3)   # root exceeds Float64 range
    @test IMU._iroot_nonneg(big(2)^999, 999) == iroot(big(2)^999, 999)  # huge x, tiny root
    @test IMU._iroot_nonneg((big(7)^200), 200) == 7                     # exact wide power
end

@testset "machine-int paths don't allocate" begin
    @test @allocated(iroot(10^15, 2)) == 0
    @test @allocated(iroot(10^15, 3)) == 0
    @test @allocated(iroot(Int128(2)^100, 5)) == 0
    @test @allocated(rootrem(10^15, 3)) == 0
    @test @allocated(ispower(10^15)) == 0
    @test @allocated(find_exponent(Int64(3)^39)) == 0
end

end
