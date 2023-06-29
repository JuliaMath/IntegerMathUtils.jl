# IntegerMathUtils.jl

This library adds several functions useful for doing math on integers. Most of these are GMP wrappers that may have faster implimentations for smaller integer types.

**Functions**

* `iroot(x::Integer, n::integer)` the integer nth root of `x`. Specifically, this is the largest integer `a` such that `a^n <= x`. Note that `n` must fit into an `Int64` (for GMP compatability).
* `ispower(x::Integer)` return if there are integer `base` and `exponent>1` values such that `base^exponent = x`.
* `find_exponent(x::Integer)` returns the largest possible integer `exponent` such that `base^exponent = x` for some `base`. Returns `1` for `x ∈ [0,1]`.
* `is_probably_prime(x::Integer; reps=25)` returns if `x` is prime. Will be incorrect less than `4^-reps` of the time.
* `kronecker(a::Integer, n::Integer)` Computes the [Kronecker_symbol](https://en.wikipedia.org/wiki/Kronecker_symbol) which is a generalization of the legendre and jacobi symbols.
