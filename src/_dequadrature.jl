
export deformulaMinusOneToOne, deformulaZeroToInf, deformulaMinusInfToInf, deint

"""
    Formula{T,P1,P2}

A variable transform for the double exponential formula. `range` gives the
truncation bounds on the transformed domain, `phi` maps a transformed node `t`
to the original domain, and `phidash` is its derivative `dx/dt`.

`phi` and `phidash` are stored with their concrete types so that node
evaluation is statically dispatched.
"""
struct Formula{T <: Real, P1, P2}
    range::Tuple{T,T}
    phi::P1
    phidash::P2
end

"""
    deformula_zero_to_inf(::Type{T}=Float64) where {T<:Real}

Constructor for the semi-infinite range (0, Inf) DE formula specialized to
numeric type `T`.
"""
deformula_zero_to_inf(::Type{T}=Float64) where {T<:Real} = Formula(
    (-T(6.8), T(6.8)),
    t -> exp(T(pi) * sinh(t) / T(2)),
    t -> T(pi) * cosh(t) * exp(T(pi) * sinh(t) / T(2)) / T(2)
)

const deformulaZeroToInf = deformula_zero_to_inf(Float64)

"""
    deformula_minus_one_to_one(::Type{T}=Float64) where {T<:Real}

Constructor for the finite range [-1, 1] DE formula specialized to numeric type `T`.
"""
deformula_minus_one_to_one(::Type{T}=Float64) where {T<:Real} = Formula(
    (-T(3.0), T(3.0)),
    t -> tanh(T(pi) * sinh(t) / T(2)),
    t -> begin
        pisinh2 = T(pi) * sinh(t) / T(2)
        sech = one(T) / cosh(pisinh2)
        T(pi) * cosh(t) * sech * sech / T(2)
    end
)

const deformulaMinusOneToOne = deformula_minus_one_to_one(Float64)

"""
    deformula_minus_inf_to_inf(::Type{T}=Float64) where {T<:Real}

Constructor for the doubly infinite range (-Inf, Inf) DE formula specialized to
numeric type `T`.

The truncation bound 6.8 is chosen so that `phidash` stays finite in `Float64`;
beyond it `cosh(pi*sinh(t)/2)` overflows.
"""
deformula_minus_inf_to_inf(::Type{T}=Float64) where {T<:Real} = Formula(
    (-T(6.8), T(6.8)),
    t -> sinh(T(pi) * sinh(t) / T(2)),
    t -> T(pi) * cosh(t) * cosh(T(pi) * sinh(t) / T(2)) / T(2)
)

const deformulaMinusInfToInf = deformula_minus_inf_to_inf(Float64)

"""
    _calcWeight!(data, t, f, phi, phidash; dropzero = zero(T))

Evaluate the integrand at the transformed node `t` and append
`(t, x, w)` to `data`.

Nodes whose weight is `NaN` are discarded. This is what lets integrable endpoint
singularities (e.g. `x^-0.5`) be handled at all: at the extreme nodes the
transform produces `phidash = Inf` while the integrand underflows to zero, and
the resulting `Inf * 0` is not a usable contribution.

Nodes whose weight magnitude does not exceed `dropzero` are discarded as well.
That is an optional node-count optimization, not a correctness requirement.
"""
function _calcWeight!(data::Vector{Tuple{T,T,T}}, t::T, f, phi, phidash; dropzero::T = zero(T)) where {T <: Real}
    xtmp::T = phi(t)
    wtmp::T = phidash(t) * f(xtmp)
    if !isnan(wtmp) && abs(wtmp) > dropzero
        if !isfinite(wtmp)
            error("Error: weight is not finite (w = $wtmp at t = $t); the integrand may diverge on the transformed domain")
        end
        push!(data, (t, xtmp, wtmp))
    end
end

"""
    _deint(f, formula; reltol::T = 1.0e-8, abstol::T = zero(T), dropzero::T = zero(T), d = 8, maxiter = 12)

Compute the numerical integration for f with double exponential formula.

    int f(x) dx = h * sum w_i

Parameters:
- f: integrand function
- formula: double expoential formula.
- reltol: tolerance for relative errors
- abstol: tolerance for absolute errors. Zero by default so that `reltol` alone
  decides convergence; a non-zero value is an early exit that caps the attainable
  relative accuracy at `abstol / |integral|`.
- dropzero: nodes whose weight magnitude does not exceed this value are discarded.
  Zero by default. Since this threshold is absolute, a non-zero value truncates
  integrands whose overall magnitude is comparable to it.
- d: the initial number of divides
- maxiter: the maximum number of iterations to increase the number of divides twice.
Return value (tuple):
- s: the value of integration
- t: a sequence for divides on the transformed domain
- x: a sequence for divides; x_i
- w: a sequence of weights (unscaled); w_i
- h: a scale parameter for w_i
"""

function _deint(f, formula::Formula{T};
        reltol::T = T(1.0e-8), abstol::T = zero(T), dropzero::T = zero(T),
        d = 8, maxiter = 12) where {T <: Real}
    lower::T = formula.range[1]
    upper::T = formula.range[2]
    h::T = (upper - lower) / d

    # Only the initial nodes are reserved up front; `push!` grows the buffer as
    # refinement actually happens. Sizing for the maxiter worst case here would
    # allocate megabytes for integrands that converge in a handful of nodes.
    data = Vector{Tuple{T,T,T}}()
    sizehint!(data, d + 1)

    for t = LinRange(lower, upper, d+1)
        _calcWeight!(data, t, f, formula.phi, formula.phidash, dropzero=dropzero)
    end

    # Compute sum directly without intermediate array
    s::T = zero(T)
    for item in data
        s += item[3]
    end
    s *= h

    iter = 1
    prev::T = s
    while true
        iter += 1
        if iter > maxiter
            @warn "Max iterations reached in _deint; returning last estimate" maxiter=maxiter s=s h=h d=d
            break
        end
        h /= 2
        for t = LinRange(lower+h, upper-h, d)
            _calcWeight!(data, t, f, formula.phi, formula.phidash, dropzero=dropzero)
        end
        d *= 2

        # Compute sum directly
        s = zero(T)
        for item in data
            s += item[3]
        end
        s *= h

        aerror = abs(s - prev)
        # The relative test must stay relative: flooring the denominator at one
        # would silently degrade it into an absolute test for integrals whose
        # magnitude is below 1.
        denom = max(abs(s), abs(prev), floatmin(T))
        rerror = aerror / denom
        if (aerror <= abstol) || (rerror <= reltol)
            break
        end
        prev = s
    end
    sort!(data, by=first)

    # Pre-allocate output vectors
    n = length(data)
    t_vec = Vector{T}(undef, n)
    x_vec = Vector{T}(undef, n)
    w_vec = Vector{T}(undef, n)
    @inbounds for i in 1:n
        t_vec[i] = data[i][1]
        x_vec[i] = data[i][2]
        w_vec[i] = data[i][3]
    end

    (s=s, t=t_vec, x=x_vec, w=w_vec, h=h)
end

"""
    deint(f, lower, upper; reltol, abstol, dropzero, d, maxiter)

Apply the double exponential (tanh-sinh) procedure to `f` on [`lower`, `upper`].

Both `lower` and `upper` may be infinite; the finite, semi-infinite (on either
side) and doubly infinite cases are dispatched onto the corresponding DE
mapping.

This function returns not only the integral value, but also the nodes and
weights used in the computation. It is intended for research and analysis
purposes rather than as a black-box integrator.

Returns a named tuple `(s, t, x, w, h)` with the integral estimate, nodes in the
transformed domain, mapped points, weights, and scale factor:
- `s`: integral estimate
- `t`: nodes in the transformed domain, in ascending order
- `x`: mapped nodes in the original domain, in the order induced by `t`
- `w`: unscaled weights corresponding to `x`
- `h`: scale factor applied to weights, so that `h * sum(w) == s`
"""
function deint(f, lower::T, upper::T;
    reltol::T=T(1.0e-8), abstol::T=zero(T), dropzero::T=zero(T),
    d=8, maxiter=12) where {T<:AbstractFloat}

    if isnan(lower) || isnan(upper)
        throw(DomainError((lower, upper), "integration bounds must not be NaN"))
    end
    if isinf(lower) && lower > zero(T)
        throw(DomainError(lower, "lower bound must not be +Inf"))
    end
    if isinf(upper) && upper < zero(T)
        throw(DomainError(upper, "upper bound must not be -Inf"))
    end

    opts = (reltol=reltol, abstol=abstol, dropzero=dropzero, d=d, maxiter=maxiter)

    if isinf(lower)
        if isinf(upper)
            # (-Inf, Inf)
            return _deint(f, deformula_minus_inf_to_inf(T); opts...)
        end
        # (-Inf, upper]: substitute x = upper - u with u in (0, Inf)
        result = _deint(deformula_zero_to_inf(T); opts...) do u
            f(upper - u)
        end
        x = upper .- result.x
        return (s=result.s, t=result.t, x=x, w=result.w, h=result.h)
    elseif isinf(upper)
        # [lower, Inf)
        if lower == zero(T)
            return _deint(f, deformula_zero_to_inf(T); opts...)
        end
        result = _deint(deformula_zero_to_inf(T); opts...) do x
            f(x + lower)
        end
        x = result.x .+ lower
        return (s=result.s, t=result.t, x=x, w=result.w, h=result.h)
    else
        # Finite interval
        if lower == -one(T) && upper == one(T)
            return _deint(f, deformula_minus_one_to_one(T); opts...)
        end
        d_half = (upper - lower) / T(2)
        result = _deint(deformula_minus_one_to_one(T); opts...) do x
            f(d_half * (x + one(T)) + lower) * d_half
        end
        x = @. d_half * (result.x + one(T)) + lower
        return (s=result.s, t=result.t, x=x, w=result.w, h=result.h)
    end
end

# Promote mixed or integer bounds, so that e.g. `deint(f, 0, Inf)` works.
function deint(f, lower::Real, upper::Real; kwargs...)
    deint(f, promote(float(lower), float(upper))...; kwargs...)
end
