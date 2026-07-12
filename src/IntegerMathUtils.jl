module IntegerMathUtils
export iroot, ispower, rootrem, find_exponent, is_probably_prime, kronecker

using Base: uabs, BitUnsigned64

@static if isdefined(Base, :top_set_bit)
    using Base: top_set_bit
else
    top_set_bit(x::Base.BitInteger) = 8sizeof(x) - leading_zeros(x)
    top_set_bit(x::Integer) = ndigits(x, base=2)
end

function iroot(x::BigInt, n::Integer)
    n <= 0 && throw(DomainError(n, "Exponent must be > 0"))
    x < 0 && iseven(n) && throw(DomainError(x, "Even roots require x >= 0"))
    ans = BigInt()
    @ccall :libgmp.__gmpz_root(ans::Ref{BigInt}, x::Ref{BigInt}, n::Culong)::Cint
    ans
end

function rootrem(x::BigInt, n::Integer)
    n <= 0 && throw(DomainError(n, "Exponent must be > 0"))
    x < 0 && iseven(n) && throw(DomainError(x, "Even roots require x >= 0"))
    root = BigInt()
    rem  = BigInt()
    @ccall :libgmp.__gmpz_rootrem(root::Ref{BigInt}, rem::Ref{BigInt}, x::Ref{BigInt}, n::Culong)::Cvoid
    return (root, rem)
end

function ispower(x::BigInt)
    return 0 != @ccall :libgmp.__gmpz_perfect_power_p(x::Ref{BigInt})::Cint
end

# a^n in T via power-by-squaring, detecting overflow. Returns (value, overflowed)
function pow_with_overflow(a::T, n::Integer) where T<:Base.BitInteger
    p = one(T)
    while n > 0
        if isodd(n)
            p, o = Base.Checked.mul_with_overflow(p, a)
            o && return (p, true)
        end
        n >>= 1
        n > 0 || break
        a, o = Base.Checked.mul_with_overflow(a, a)
        o && return (p, true)
    end
    return (p, false)
end
function _pow_leq(a::T, n::Integer, x::T) where T<:Base.BitInteger
    p, o = pow_with_overflow(a, n)
    return !o && p <= x
end
_pow_leq(a::T, n::Integer, x::T) where T<:Integer = a^n <= x

# Float64-accurate estimate of x^(1/n) for x <= typemax(UInt128).
# Because the root must be small, this is within +-1 of the true root when n>2.
function _iroot_approx(x::T, n::Integer) where T<:Base.BitInteger
    trunc(T, exp(log(Float64(x)) / n))
end

# A Float64-accurate estimate of x^(1/n) for x > typemax(UInt128)
# Split x = hi*2^sh with hi its top 53 bits, estimate log2(root) = (log2(hi)+sh)/n
# and place a 53-bit mantissa at the right exponent.
function _iroot_float(x::T, n::Integer) where T<:Integer
    b = top_set_bit(x)
    sh = b - 53
    L = (log2(Float64(x >> sh)) + sh) / n
    e = floor(Int, L)
    m53 = trunc(T, exp2(L - e) * exp2(53))
    return e >= 53 ? (m53 << (e - 53)) : (m53 >> (53 - e))
end

# One integer Newton step
@inline function _iroot_newton(x::T, n::Integer, a::T) where T<:Integer
    return ((n - 1) * a + x ÷ a^(n - 1)) ÷ n
end

function _iroot_approx(x::T, n::Integer) where T<:Integer
    x <= typemax(UInt128) && return T(_iroot_approx(UInt128(x), n))
    # one unconditional step lifts any seed to >= floor(root)-1 (by weighted AM-GM)
    a = _iroot_newton(x, n, _iroot_float(x, n))
    while true # after which the iteration descends until return
        b = _iroot_newton(x, n, a)
        b >= a && return a
        a = b
    end
end

# floor(x^(1/n)) for x >= 0; n >= 2 (n == 2 only for non-bit Integer types,
# machine ints take the Base.isqrt path in iroot)
function _iroot_nonneg(x::T, n::Integer) where T<:Integer
    x <= one(T) && return x
    n >= top_set_bit(x) && return one(T)   # 2^n > x
    a = _iroot_approx(x, n)
    if _pow_leq(a, n, x)
        while _pow_leq(a + one(T), n, x)
            a += one(T)
        end
    else
        a -= one(T)
        while !_pow_leq(a, n, x)
            a -= one(T)
        end
    end
    return a
end

# Base.isqrt is already optimal for machine ints; its generic Integer
# fallback is Float64-based though, which is wrong for wide values.
_isqrt(x::Base.BitInteger) = Base.isqrt(x)
_isqrt(x::Integer) = _iroot_nonneg(x, 2)

function iroot(x::T, n::Integer) where T<:Integer
    n <= 0 && throw(DomainError(n, "Exponent must be > 0"))
    x < 0 && iseven(n) && throw(DomainError(x, "Even roots require x >= 0"))
    n == 1 && return x
    n == 2 && return T(_isqrt(uabs(x)))   # x >= 0 here
    r = T(_iroot_nonneg(uabs(x), n))
    return x < 0 ? -r : r
end

function rootrem(x::Integer, n::Integer)
    a = iroot(x, n)
    return (a, x - a^n)
end

# Is m (>= 2) an exact k-th power (k odd, >= 3)?
# Callers handle k == 2 via _isqrt directly
# (its root can be large enough that the float error below exceeds 0.5).
# The float root r = exp(lm/k) has provable worst-case error
# E <= 1.5*r*(|ln r|+1)*2^-52; for k >= 3 and m < 2^128 this is << 0.5
# so round(r) is the only integer candidate
function _kth_root_exact(m::T, k::Integer, logm::Float64) where T<:Base.BitInteger
    lr = logm / k
    r = exp(lr)
    tol = 2 * r * (abs(lr) + 1) * exp2(-52)
    b = round(r)
    abs(r - b) > tol && return false
    c = unsafe_trunc(T, b)
    p, o = pow_with_overflow(c, k)
    return !o && p == m
end
# Arbitrary-precision types: no reliable Float64, so verify via the native root.
function _kth_root_exact(m::T, k::Integer, logm::Float64) where T<:Integer
    a = iroot(m, k)
    return a^k == m # a^k <= m, so no overflow
end

function ispower(x::Integer)
    ux = uabs(x)
    ux <= one(ux) && return true
    neg = x < 0
    trailing_zeros(ux) == 1 && return false # 2 | x exactly once: not a power
    if !neg
        a = _isqrt(ux)
        a * a == ux && return true
    end
    lx = log(Float64(ux))
    for k in 3:2:top_set_bit(ux)
        _kth_root_exact(ux, k, lx)[1] && return true
    end
    return false
end

function find_exponent(x::Integer)
    ux = uabs(x)
    ux <= one(ux) && return 1
    neg = x < 0
    lx = log(Float64(ux))
    for e in top_set_bit(ux):-1:3
        neg && iseven(e) && continue
        _kth_root_exact(ux, e, lx)[1] && return e
    end
    if !neg
        a = _isqrt(ux)
        a * a == ux && return 2
    end
    return 1
end

function is_probably_prime(x::Integer; reps=25)
    if !(x isa BigInt)
        x = BigInt(x)
    end
    return 0 != @ccall :libgmp.__gmpz_probab_prime_p(x::Ref{BigInt}, reps::Cint)::Cint
end

function kronecker(a::BigInt, b::Clong)
    return @ccall :libgmp.__gmpz_kronecker_si(a::Ref{BigInt}, b::Clong)::Cint
end
function kronecker(a::Clong, b::BigInt)
    return @ccall :libgmp.__gmpz_si_kronecker(a::Clong, b::Ref{BigInt})::Cint
end
function kronecker(a, n)
    @assert n != -n || n == 0
    @assert a != -a || a == 0
    t = 1
    if iszero(n)
        return Int(abs(a) == 1)
    end
    if n < 0
        n = abs(n)
        if a < 0
            t = -t
        end
    end
    trail = trailing_zeros(n)
    if trail > 0
        n >>= trail
        if iseven(a)
            return 0
        elseif isodd(trail) && a&7 in (3,5)
            t = -t
        end
    end
    a = mod(a, n)
    while a != 0
        while iseven(a)
            a = a >> 1
            if n&7 in (3, 5)
                t = -t
            end
        end
        a, n = n, a
        if a&3 == n&3 == 3
            t = -t
        end
        a = mod(a, n)
    end
    return n == 1 ? t : 0
end

end
