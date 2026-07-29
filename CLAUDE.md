# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

DEQuadrature.jl is a small, dependency-free Julia package implementing double-exponential (tanh–sinh)
quadrature on finite, semi-infinite and doubly infinite intervals. It returns not just the integral
value but the full node/weight set, so it is meant as an analysis tool, not a black-box integrator.

The package was renamed from `Deformula` to `DEQuadrature`; the working directory is still
`Deformula.jl` and some testset names still say `Deformula*`.

## Release and registry

The package is **not** in the General registry. It is registered in the org's own registry,
`github.com/JuliaReliab/Registry` (entry `D/DEQuadrature`), which is managed with LocalRegistry.jl.

To release a new version:

1. Bump `version` in `Project.toml` and add a `CHANGELOG.md` section.
2. Merge to `master`. Note `master` requires one approving review, and GitHub does not allow
   self-approval, so a solo merge needs `gh pr merge <n> --merge --admin`.
3. Tag and release manually — `TagBot` only works for General registrations:
   `git tag -a vX.Y.Z -m "DEQuadrature X.Y.Z" && git push origin vX.Y.Z`, then `gh release create`.
4. Register. LocalRegistry requires a **clean working tree** and reads the tree hash of `HEAD`, which
   must already be pushed:

```julia
using LocalRegistry   # not a dep of this project; load it from another environment
register("/Users/okamu/Documents/Deformula.jl", registry = "Registry")
```

Verify in a temp environment (`Pkg.activate(mktempdir())`) so the `path = "."` dev entry cannot mask
a broken registration.

If DEQuadrature is ever submitted to General, AutoMerge will flag the name against the existing
`Quadrature` package: Damerau–Levenshtein distance must be ≥ 3 and these are 2 apart. That needs the
`Override AutoMerge: name similarity is okay.` label or a rename.

## Commands

```julia
# Run the whole test suite
julia --project=. -e 'using Pkg; Pkg.test()'
```

Julia has no per-testset filter. To iterate on one testset, temporarily comment out the others in
`test/runtests.jl`, or copy the `@testset` body into a REPL after
`import DEQuadrature: _deint, deint, deformulaZeroToInf, deformulaMinusOneToOne, deformulaMinusInfToInf`.

CI (`.github/workflows/ci.yml`) tests Julia 1.10, `1`, and nightly on Linux/macOS/Windows with
coverage → Codecov. `Project.toml` declares `julia = "1.6"` compat, which is looser than what CI
exercises.

## Architecture

Two files only: `src/DEQuadrature.jl` is a shell that `include`s `src/_dequadrature.jl`, which holds
everything.

The design is a two-layer split:

- **`Formula{T,P1,P2}`** — a variable transform bundling `range::Tuple{T,T}` (truncation bounds in the
  transformed `t` domain) with `phi` (t → x) and `phidash` (dx/dt). `phi`/`phidash` are type
  parameters on purpose: leaving them untyped made every node evaluation a dynamic dispatch and cost
  roughly 2× throughput. Three built-in mappings: `deformula_zero_to_inf(T)` for (0, ∞) with range
  ±6.8, `deformula_minus_one_to_one(T)` for (−1, 1) with range ±3.0, and
  `deformula_minus_inf_to_inf(T)` for (−∞, ∞) with range ±6.8. The exported `deformulaZeroToInf` /
  `deformulaMinusOneToOne` / `deformulaMinusInfToInf` are the `Float64` instances.
- **`_deint(f, formula)`** — the trapezoid-with-halving engine. It evaluates `d+1` nodes, then each
  iteration halves `h` and adds only the *new* midpoints to a shared `data::Vector{Tuple{t,x,w}}`,
  doubling `d`. Because nodes accumulate rather than being recomputed, `data` is sorted by `t` only
  at the end and unpacked into the returned `t`/`x`/`w` vectors. Convergence is
  `aerror <= abstol || rerror <= reltol`; hitting `maxiter` emits a `@warn` and returns the last
  estimate.
- **`deint(f, lower, upper)`** — dispatches an arbitrary interval onto one of the three formulas:
  `(-∞,∞)` direct, `(-∞,b]` by the substitution `x = b - u`, `[0,∞)` direct, `[lower,∞)` by shifting
  the integrand, `[-1,1]` direct, and any other finite `[a,b]` by affine rescaling with the Jacobian
  `d_half` folded into the integrand. `+Inf` as a lower bound, `-Inf` as an upper bound and `NaN`
  throw `DomainError`. A second method promotes mixed/integer bounds so `deint(f, 0, Inf)` works.
  In the substituted, shifted and rescaled branches only `x` is mapped back — `w` and `h` stay in the
  reference domain, and the identity `h * sum(w) ≈ s` still holds. `x` follows the order induced by
  `t`, so it is *descending* for `(-∞,b]` and for reversed finite intervals.

## Invariants worth not breaking

These each encode a bug that was fixed in 0.3.0, so a "cleanup" that reverts one silently returns
wrong numbers.

- **`_calcWeight!` must keep dropping NaN weights.** At the extreme nodes the transform produces
  `phidash = Inf` while the integrand underflows to 0, and `Inf * 0` is not a usable contribution.
  This is what makes integrable endpoint singularities (e.g. `x^-0.5`) work at all.
- **`abstol` and `dropzero` are separate, and both default to `zero(T)`.** `dropzero` is an
  *absolute* node-dropping threshold; when it defaulted to `eps(T)` (and was spelled `abstol`),
  integrands of comparable magnitude lost most of their nodes — `∫₀^∞ 1e-300 e^{-x} dx` returned
  exactly `0.0`. Likewise a non-zero `abstol` is an early exit that caps attainable relative accuracy
  at `abstol / |integral|`.
- **The relative test's denominator must not be floored at one.** `max(abs(s), abs(prev),
  floatmin(T))` is deliberate; flooring at `one(T)` degrades the relative check into an absolute one
  for integrals of magnitude below 1.
- **`sizehint!` reserves only `d+1`.** Sizing for the `maxiter` worst case allocated 6 MB at
  `maxiter=16` for a 129-node result, with no accuracy benefit.

Known limitation, not a bug to chase: the error estimate is the difference between successive
refinements and can be optimistic. `∫₀¹ x^-0.5 dx` converges to ~1.4e-7 relative error against
`reltol=1e-8`. Fixing it needs a different error estimator.

## Calling convention gotcha

`f` is the *first* positional argument of both `_deint` and `deint`, so `do`-block form reads
`_deint(deformulaZeroToInf) do x ... end` and `deint(0.0, Inf64) do x ... end` — the formula/bounds
appear alone in the call. Tests use this form throughout.

## BigFloat support

There is a single generic `deint(f, ::T, ::T) where {T<:AbstractFloat}` (the `Float64`-specialized
duplicate was removed in 0.3.0, along with its conflicting `reltol` default). All numeric literals go
through `T(...)`, so `Formula`/`_deint` work at any precision. Use `setprecision(n) do ... end` and
pass `BigFloat` bounds plus a tightened `reltol`. Keep new code type-generic in `T`; the only
intentionally `Float64`-bound pieces are the three exported consts.
