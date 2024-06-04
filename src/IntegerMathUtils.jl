module IntegerMathUtils
export iroot, ispower, rootrem, find_exponent, is_probably_prime, kronecker

iroot(x::Integer, n::Integer) = iroot(x, Cint(n))

# TODO: Add more efficient implimentation
iroot(x::T, n::Cint) where {T<:Integer} = T(iroot(big(x), Cint(n)))

function iroot(x::BigInt, n::Cint)
    n <= 0 && throw(DomainError(n, "Exponent must be > 0"))
    x <= 0 && iseven(x) && throw(DomainError(n, "This is a math no-no"))
    ans = BigInt()
    @ccall :libgmp.__gmpz_root(ans::Ref{BigInt}, x::Ref{BigInt}, n::Cint)::Cint
    ans
end

# TODO: Add more efficient implimentation for smaller types
function rootrem(x::T, n::Integer) where {T<:Integer}
    x = big(x)
    n = Cint(n)
    n <= 0 && throw(DomainError(n, "Exponent must be > 0"))
    x <= 0 && iseven(x) && throw(DomainError(n, "This is a math no-no"))
    root = BigInt()
    rem  = BigInt()
    @ccall :libgmp.__gmpz_rootrem(root::Ref{BigInt}, rem::Ref{BigInt}, x::Ref{BigInt}, n::Int)::Nothing
    return (root, T(rem))
end

# TODO: Add more efficient implimentation for smaller types
ispower(x::Integer) = ispower(big(x))

function ispower(x::BigInt)
    return 0 != @ccall :libgmp.__gmpz_perfect_power_p(x::Ref{BigInt})::Cint
end

# TODO: Add more efficient implimentation for smaller types
function find_exponent(x::Integer)
    x <= 1 && return 1
    for exponent in ndigits(x, base=2):-1:2
        rootrem(x, exponent)[2] == 0 && return exponent
    end
    1
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

struct Montgomery{T<:Integer}
    n::T
    nr::T
end
function Montgomery(n::T) where T
    nr = one(T)
    num_iters = trailing_zeros(8*sizeof(T))
    for i in 1:num_iters
        nr *= 2 - n * nr
    end
    return Montgomery{T}(n, nr)
end
function montgomify(M::Montgomery{T}, x::T) where T
    shift = 8*sizeof(T)
    return T((widen(x) << shift) % M.n)
end
# returns a number in the [0, 2 * n - 2] range
function montgomery_reduce(M::Montgomery{T}, x) where T
    shift = 8*sizeof(T)
    q = (x%T) * M.nr
    m = widemul(q, M.n) >> shift
    return ((x >> shift) + M.n - m) % T
end
function montgomery_mul(M::Montgomery{T}, x::T, y::T) where T
    return montgomery_reduce(M, widemul(x, y))
end
function Base.powermod(a::T, p::T, M::Montgomery{T}) where T
    a = montgomify(M, a)
    r = montgomify(M,one(T))
    while p != 0
        if p&1 != 0
            r = montgomery_mul(M, r, a)
        end
        a = montgomery_mul(M, a, a)
        p >>= 1
    end
    ans = montgomery_reduce(M, r)
    return ans < M.n ? ans : ans - M.n
end

end
