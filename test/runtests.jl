import DEQuadrature: _deint, deformulaZeroToInf, deformulaMinusOneToOne,
                     deformulaMinusInfToInf, deint
using Test

@testset "DEQuadrature.jl" begin
    result = _deint(deformulaZeroToInf) do x
        0.5*x^-0.5*exp(-x^0.5)
    end
    @test result.s ≈ 1.0
    result = _deint(deformulaMinusOneToOne) do x
        y = x + 1
        0.5*y^-0.5*exp(-y^0.5)
    end
    @test result.s ≈ 0.7568830631028768
end

@testset "Deformula2" begin
    result = _deint(deformulaZeroToInf) do x
        0.5*x^-0.5*exp(-x^0.5)
    end
    @test sum(result.w) * result.h ≈ 1.0
end

@testset "Deformula3" begin
    result = deint(0.0, Inf64) do x
        0.5*x^-0.5*exp(-x^0.5)
    end
    @test sum(result.w) * result.h ≈ 1.0
end

@testset "Deformula4" begin
    result = deint(0.0, 10.0) do x
        if 0 <= x <= 10
            1
        else
            0
        end
    end
    @test result.s ≈ 10.0
end

@testset "Deformula5" begin
    result = deint(1.0, Inf64) do x
        1.0/sqrt(2.0*pi) * exp(-(x - 1.0)^2/2.0)
    end
    @test result.s ≈ 0.5
end

# Additional tests

@testset "FiniteIntervalGeneral" begin
    # Integrate constant 1 over a non-symmetric finite interval
    a, b = -2.5, 3.7
    result = deint(a, b) do x
        1.0
    end
    @test result.s ≈ (b - a) atol=1e-8 rtol=1e-8
    @test issorted(result.t)
end

@testset "BigFloatSemiInfinite" begin
    setprecision(256) do
        result = deint(BigFloat(0), BigFloat(Inf); reltol=BigFloat(1e-20)) do x
            exp(-x)
        end
        @test result.s ≈ BigFloat(1) rtol=BigFloat(1e-18)
    end
end

@testset "BigFloatFiniteInterval" begin
    setprecision(256) do
        a = BigFloat(-2)
        b = BigFloat(3)
        result = deint(a, b; reltol=BigFloat(1e-20)) do x
            one(BigFloat)
        end
        @test result.s ≈ (b - a) rtol=BigFloat(1e-12)
        @test issorted(result.t)
    end
end

@testset "MaxIterWarning" begin
    # Force a warning by setting maxiter very small; ensure it doesn't throw
    @test_logs (:warn,) begin
        result = _deint(deformulaZeroToInf; maxiter=1) do x
            exp(-x)
        end
        @test result.s > 0
    end
end

@testset "SemiInfiniteShift" begin
    # Check [lower, Inf) with a shift works; mean=2 normal -> P(X>=2)=0.5
    result = deint(2.0, Inf64) do x
        1.0 / sqrt(2.0*pi) * exp(-((x - 2.0)^2) / 2.0)
    end
    @test result.s ≈ 0.5 atol=1e-8 rtol=1e-8
end

@testset "Deformula6" begin
    dfm(x; m, lambda) = lambda * x^(m-1)* exp(-lambda * x^m) + lambda * m * x^(m-1) * log(x) * exp(-lambda * x^m) + lambda * m * x^(m-1) * exp(-lambda * x^m) * (-lambda * x^m * log(x))
    # This integrand is the partial derivative w.r.t. m of the Weibull PDF λ m x^{m-1} e^{-λ x^m}.
    # Therefore, ∫_0^∞ dfm(x; m, λ) dx = ∂/∂m ∫_0^∞ PDF dx = ∂/∂m 1 = 0.
    # We assert the integral is approximately 0 within absolute tolerance.
    result = deint(0.0, Inf64; reltol=1e-8, abstol=1e-12, d=8, maxiter=20) do x
        dfm(x, m=2.0, lambda=1.0)
    end
    @test isfinite(result.s)
    @test result.s ≈ 0.0 atol=1e-8
end

# Regression tests

@testset "InfiniteLowerBound" begin
    # Used to return 0.0 silently: the -Inf lower bound fell through to the
    # finite-interval branch, where d_half = Inf killed every weight.
    result = deint(-Inf, 0.0) do x
        exp(x)
    end
    @test result.s ≈ 1.0 rtol=1e-8
    @test result.h * sum(result.w) ≈ result.s

    result = deint(-Inf, 2.0) do x
        exp(x)
    end
    @test result.s ≈ exp(2.0) rtol=1e-8
    @test all(result.x .<= 2.0)
end

@testset "DoublyInfinite" begin
    result = deint(-Inf, Inf) do x
        exp(-x^2)
    end
    @test result.s ≈ sqrt(pi) rtol=1e-8
    @test result.h * sum(result.w) ≈ result.s

    result = deint(-Inf, Inf) do x
        1.0 / (1.0 + x^2)
    end
    @test result.s ≈ pi rtol=1e-8

    result = _deint(deformulaMinusInfToInf) do x
        exp(-x^2)
    end
    @test result.s ≈ sqrt(pi) rtol=1e-8
    @test issorted(result.t)
end

@testset "InvalidBounds" begin
    @test_throws DomainError deint(x -> 1.0, Inf, Inf)
    @test_throws DomainError deint(x -> 1.0, 0.0, -Inf)
    @test_throws DomainError deint(x -> 1.0, NaN, 1.0)
    @test_throws DomainError deint(x -> 1.0, 0.0, NaN)
end

@testset "SmallMagnitudeRelativeAccuracy" begin
    # The relative test used to floor its denominator at 1, and the default
    # node-dropping threshold was an absolute eps(T). Both turned the requested
    # relative tolerance into an absolute one for small-valued integrals; the
    # 1e-300 case returned exactly 0.0.
    for scale in (1e-13, 1e-100, 1e-300)
        result = deint(0.0, Inf64) do x
            scale * exp(-x)
        end
        @test result.s ≈ scale rtol=1e-8
    end
end

@testset "LargeMaxIterDoesNotOverAllocate" begin
    # sizehint! used to reserve the maxiter worst case up front, so raising
    # maxiter inflated memory by orders of magnitude without adding a node.
    f = x -> 0.5*x^-0.5*exp(-x^0.5)
    r12 = _deint(f, deformulaZeroToInf, maxiter=12)
    r20 = _deint(f, deformulaZeroToInf, maxiter=20)
    @test length(r12.x) == length(r20.x)
    a12 = @allocated _deint(f, deformulaZeroToInf, maxiter=12)
    a20 = @allocated _deint(f, deformulaZeroToInf, maxiter=20)
    @test a20 < 2 * a12
end

@testset "DropzeroIsSeparateFromAbstol" begin
    # abstol must not silently discard nodes any more.
    exact = deint(x -> exp(-x), 0.0, Inf64)
    loose = deint(x -> exp(-x), 0.0, Inf64, abstol=1e-3)
    @test loose.s ≈ 1.0 rtol=1e-5
    # dropzero is the knob that trades nodes for accuracy.
    dropped = deint(x -> exp(-x), 0.0, Inf64, dropzero=1e-3)
    @test isapprox(dropped.s, 1.0; rtol=1e-3)
    @test all(abs.(dropped.w) .> 1e-3)
end

@testset "PromotedBounds" begin
    @test deint(x -> exp(-x), 0, Inf).s ≈ 1.0 rtol=1e-8
    @test deint(x -> 2x, 0, 1).s ≈ 1.0 rtol=1e-8
    @test deint(x -> 1.0, 0, 10).s ≈ 10.0 rtol=1e-8
end

@testset "BigFloatDoublyInfinite" begin
    setprecision(128) do
        result = deint(BigFloat(-Inf), BigFloat(Inf); reltol=BigFloat(1e-20)) do x
            exp(-x^2)
        end
        @test result.s ≈ sqrt(BigFloat(pi)) rtol=BigFloat(1e-18)
    end
end

@testset "FormulaIsConcretelyTyped" begin
    # Untyped phi/phidash fields made every node evaluation a dynamic dispatch.
    F = deformulaZeroToInf
    @test isconcretetype(fieldtype(typeof(F), :phi))
    @test isconcretetype(fieldtype(typeof(F), :phidash))
end
