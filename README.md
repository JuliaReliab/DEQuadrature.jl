# DEQuadrature.jl

[![CI](https://github.com/JuliaReliab/DEQuadrature.jl/actions/workflows/ci.yml/badge.svg?branch=master)](https://github.com/JuliaReliab/DEQuadrature.jl/actions/workflows/ci.yml)
[![codecov](https://codecov.io/gh/JuliaReliab/DEQuadrature.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JuliaReliab/DEQuadrature.jl)

Fast, robust numerical integration on finite and semi-infinite intervals using the Double Exponential (tanh–sinh) quadrature.

DEQuadrature.jl provides a simple, allocation‑light API to integrate 1‑D functions with good accuracy for smooth and rapidly decaying integrands.

## Installation

Once registered in the General registry:

```julia
using Pkg
Pkg.add("DEQuadrature")
```

Until then, you can install from GitHub:

```julia
using Pkg
Pkg.add(url = "https://github.com/JuliaReliab/DEQuadrature.jl.git")
```

## Quick start

```julia
using DEQuadrature

# Integrate e^{-x} on [0, ∞)
res = deint(x -> exp(-x), 0.0, Inf64)
@show res.s           # integral value
@show res.h * sum(res.w)  # scaled weights (sanity check)

# Integrate on a finite interval [a,b]
a, b = -2.5, 3.7
res2 = deint(x -> 1.0, a, b)
@show res2.s  # ≈ b - a

# Semi-infinite on the left, and doubly infinite
@show deint(x -> exp(x), -Inf, 0.0).s      # ≈ 1
@show deint(x -> exp(-x^2), -Inf, Inf).s   # ≈ √π

# Integer or mixed bounds are promoted
@show deint(x -> exp(-x), 0, Inf).s        # ≈ 1
```

## API

```julia
deint(f, lower::T, upper::T;
      reltol::T = T(1.0e-8),
      abstol::T = zero(T),
      dropzero::T = zero(T),
      d::Int = 8,
      maxiter::Int = 12) where {T<:AbstractFloat}
```

`lower` and `upper` may each be infinite. Finite, semi-infinite (on either side)
and doubly infinite intervals are dispatched onto the corresponding DE mapping;
`+Inf` as a lower bound, `-Inf` as an upper bound, and `NaN` throw a `DomainError`.
A promoting method accepts mixed or integer bounds.

Computes

$$
\int_{lower}^{upper} f(x)\,dx \;\approx\; h\,\sum_i w_i
$$

Return value is a named tuple `(s, t, x, w, h)`:

- `s`: integral estimate
- `t`: nodes in the transformed domain
- `x`: nodes in the original domain
- `w`: unscaled weights
- `h`: scale for weights (so that `h*sum(w) ≈ s`)

Notes:

- Convergence stops when either the absolute or the relative change between two
  successive refinements falls within tolerance. `abstol` defaults to zero, so by
  default `reltol` alone decides — a non-zero `abstol` is an early exit that caps
  the attainable relative accuracy at `abstol / |integral|`.
- `dropzero` discards nodes whose weight magnitude does not exceed it. The
  threshold is absolute, so a non-zero value truncates integrands whose overall
  magnitude is comparable to it. It defaults to zero for that reason; raise it
  only to trade accuracy for fewer nodes.
- On reaching `maxiter` before convergence, a warning is emitted and the last
  estimate is returned.
- The error estimate is the difference between successive refinements, which can
  be optimistic. For integrands with strong endpoint singularities (e.g. `x^-0.5`
  on `[0,1]`) the achieved relative error may be an order of magnitude larger
  than `reltol`.

### Arbitrary precision (BigFloat)

The single `deint` method above is generic in `T`, so arbitrary precision works
by passing `BigFloat` bounds.

Usage with `BigFloat`:

```julia
using DEQuadrature

setprecision(256) do
      # ∫_0^∞ e^{-x} dx = 1
      res = deint(x -> exp(-x), BigFloat(0), BigFloat(Inf); reltol=BigFloat(1e-20))
      @show res.s  # ≈ 1 with high precision

      # Finite interval [a,b]
      a = BigFloat(-2); b = BigFloat(3)
      res2 = deint(x -> one(BigFloat), a, b; reltol=BigFloat(1e-18))
      @show res2.s  # ≈ b - a
end
```

Tips:
- Set `setprecision` to suit your needs and tighten `reltol` accordingly.
- For very tight tolerances or slowly decaying integrands, consider increasing `maxiter` and/or the initial divisions `d`.

### Tuning for performance

The default parameters `(d=8, maxiter=12)` balance speed and accuracy for most smooth, well-behaved integrands. Adjust based on your needs:

**maxiter** (default: 12)
- Controls the refinement depth. Each iteration doubles the number of divisions.
- Sufficient for `reltol ~1e-8` with `Float64`. For tighter tolerances, increase modestly (e.g., to 14–16).
- Raising it costs nothing when the integrand converges early: node storage grows with the nodes actually evaluated, not with `maxiter`.

**d** (default: 8)
- Initial node count per interval. Larger `d` front-loads computation but may need fewer refinement steps.
- Try `d=10–16` for smooth integrands; `d=6–8` for rougher ones or memory-constrained settings.

**reltol** (default: `1e-8`)
- Convergence is based on the relative change between successive refinements.
- Tighten this *first* for higher accuracy; increasing `maxiter` is a fallback when tolerance targets are nearly met but time/memory is less critical.

Advanced users can also construct the DE mapping for an arbitrary `T` via:

```julia
deformula_zero_to_inf(T)        # semi-infinite (0, ∞)
deformula_minus_one_to_one(T)   # finite interval (-1, 1)
deformula_minus_inf_to_inf(T)   # doubly infinite (-∞, ∞)
```
and call the internal `_deint(f, formula::Formula; ...)` directly if needed.

### Exported symbols

- `deint` — high‑level integration API
- `deformulaZeroToInf`, `deformulaMinusOneToOne`, `deformulaMinusInfToInf` — prebuilt `Float64` DE mappings (advanced use)
      - For arbitrary precision: use `deformula_zero_to_inf(T)` / `deformula_minus_one_to_one(T)` / `deformula_minus_inf_to_inf(T)`

## Compatibility

- Julia 1.x (see `Project.toml` for details)

## Development

- Run tests:

```julia
using Pkg
Pkg.activate(".")
Pkg.test()
```

Contributions are welcome. Please open an issue or pull request.

## License

This project is licensed under the terms of the LICENSE file included in this repository.
