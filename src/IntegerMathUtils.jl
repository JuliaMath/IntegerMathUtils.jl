module IntegerMathUtils
export iroot, ispower, rootrem, find_exponent, is_probably_prime

iroot(x::Integer, n::Integer) = iroot(x, Cint(n))

# TODO: Add more efficient implimentation
iroot(x::T, n::Cint) where {T<:Integer} = T(iroot(big(x), Cint(n)))

function iroot(x::BigInt, n::Cint)
    n <= 0 && throw(DomainError(n, "Exponent must be > 0"))
    x <= 0 && iseven(x) && throw(DomainError(n, "This is a math no-no"))
    ans = BigInt()
    ccall((:__gmpz_root, :libgmp), Cint, (Ref{BigInt}, Ref{BigInt}, Cint), ans, x, n)
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
    ccall((:__gmpz_rootrem, :libgmp), Nothing,(Ref{BigInt}, Ref{BigInt}, Ref{BigInt}, Int), root, rem, x, n)
    return (root, T(rem))
end

# TODO: Add more efficient implimentation for smaller types
ispower(x::Integer) = ispower(big(x))

function ispower(x::BigInt)
    return 0 != ccall((:__gmpz_perfect_power_p, :libgmp), Cint, (Ref{BigInt},), x)
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
    return ccall((:__gmpz_probab_prime_p, :libgmp), Cint, (Ref{BigInt}, Cint), x, reps) != 0
end

end
