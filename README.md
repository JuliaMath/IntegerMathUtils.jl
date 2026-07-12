# IntegerMathUtils.jl

This library adds several functions useful for doing math on integers. Machine integers (and any other non-`BigInt` `Integer` type) use native Julia algorithms; `BigInt` inputs use GMP.

**Functions**

* `iroot(x::Integer, n::Integer)` the integer nth root of `x`, truncated toward zero: for `x >= 0` the largest integer `a` such that `a^n <= x`, and `iroot(-x, n) == -iroot(x, n)` for odd `n` (e.g. `iroot(-9, 3) == -2`). Throws a `DomainError` for `n <= 0` or an even root of a negative number.
* `rootrem(x::Integer, n::Integer)` returns `(iroot(x, n), x - iroot(x, n)^n)`.
* `ispower(x::Integer)` returns whether there are integer `base` and `exponent > 1` values such that `base^exponent == x`. For negative `x` only odd exponents count, so `ispower(-8) == true` but `ispower(-4) == false`. `0` and `±1` are perfect powers.
* `find_exponent(x::Integer)` returns the largest possible integer `exponent` such that `base^exponent == x` for some `base` (only odd exponents for negative `x`). Returns `1` for `x ∈ [-1, 1]`.
* `is_probably_prime(x::Integer; reps=25)` returns if `x` is prime. Will be incorrect less than `4^-reps` of the time.
* `kronecker(a::Integer, n::Integer)` Computes the [Kronecker symbol](https://en.wikipedia.org/wiki/Kronecker_symbol) which is a generalization of the legendre and jacobi symbols.
