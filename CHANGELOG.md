# v0.3.0

Bug fixes (all of these previously produced wrong results without any warning):

- Support infinite lower bounds. `deint(f, -Inf, b)` and `deint(f, -Inf, Inf)`
  used to fall through to the finite-interval branch and return `0.0`. Adds the
  doubly infinite mapping `deformula_minus_inf_to_inf` / `deformulaMinusInfToInf`.
- Keep the relative convergence test relative. The denominator was floored at
  one, which degraded it into an absolute test for integrals of magnitude below
  one.
- Change the default `abstol` from `eps(T)` to `zero(T)` and split the
  node-dropping threshold out of `abstol` into a new `dropzero` keyword, also
  defaulting to `zero(T)`. Because both were absolute thresholds, small-valued
  integrands lost most of their nodes; `∫_0^∞ 1e-300 e^{-x} dx` returned exactly
  `0.0`, and `1e-13 e^{-x}` was off by 0.01%.
- Reject `+Inf` as a lower bound, `-Inf` as an upper bound, and `NaN` bounds with
  a `DomainError` instead of returning a meaningless value.

Performance:

- `Formula` now stores `phi` / `phidash` with concrete types, so node evaluation
  is statically dispatched.
- `sizehint!` no longer reserves the `maxiter` worst case up front. Raising
  `maxiter` used to inflate allocation by orders of magnitude without adding a
  single node (`maxiter=16` allocated 6 MB for a 129-node result).
- Together these cut the reference `_deint` benchmark from 832 ms to 7.6 ms per
  200 runs, while allocation became independent of `maxiter`.

Other:

- Merge the `Float64`-specialized `deint` into the generic method, resolving the
  conflicting `reltol` defaults (`1e-8` / `1e-9`) between them. The public
  default is now `1e-8` for every `T`.
- Accept mixed or integer bounds, so `deint(f, 0, Inf)` works.
- Fix the non-finite weight error message, which reported `NaN` for an `Inf`.

# v0.2.0

- Change default value of maxiter from 16 to 12
- Update for BigFloat support
- Change the name of the package from Deformula to DEQuadrature

# v0.1.4

- Update compatibility to Julia 1.6 and above.
- Change ci.yml to use Julia 1.6 for testing.

# v0.1.3

- Add CI workflow
- Add CompatHelper workflow

