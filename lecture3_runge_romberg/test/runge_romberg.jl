include("../impl/runge_romberg.jl")

using Test

# ─────────────────────────────────────────────────────────────
#  Подключи свою реализацию, например:
#  include("romberg.jl")
# ─────────────────────────────────────────────────────────────

@testset "Romberg" begin

    # ── 1. Структура возвращаемых матриц ─────────────────────

    @testset "Output shape" begin
        f = x -> x^2
        U, p_eff = romberg(f, 0, 1, 10, levels=4)
        @test size(U)     == (4, 4)
        @test size(p_eff) == (4, 4)
    end

    @testset "First column is monotone convergent" begin
        # U[s,1] — трапеции с удвоением сетки, должны сходиться
        f = x -> x^2
        U, _ = romberg(f, 0, 1, 10, levels=4)
        exact = 1/3
        for s in 2:4
            @test abs(U[s, 1] - exact) < abs(U[s-1, 1] - exact)
        end
    end

    @testset "Upper triangle is zero" begin
        # U[s, l] не заполняется при l > s — должны остаться нули
        f = x -> x^2
        U, p_eff = romberg(f, 0, 1, 10, levels=4)
        for s in 1:4, l in s+1:4
            @test U[s, l]     == 0.0
            @test p_eff[s, l] == 0.0
        end
    end

    # ── 2. Точность — диагональ улучшается ───────────────────

    @testset "Diagonal improves accuracy: x^2 on [0,1]" begin
        f = x -> x^2
        exact = 1/3
        U, _ = romberg(f, 0, 1, 10, levels=4)
        # Финальный результат заметно точнее грубой трапеции
        @test U[4, 4] ≈ exact atol=1e-4
        # С большим N — высокая точность
        U2, _ = romberg(f, 0, 1, 100, levels=4)
        @test U2[4, 4] ≈ exact atol=1e-10
    end

    @testset "Diagonal improves accuracy: x^3 on [0,2]" begin
        f = x -> x^3
        exact = 4.0
        U, _ = romberg(f, 0, 2, 10, levels=4)
        @test U[4, 4] ≈ exact atol=1e-3
        U2, _ = romberg(f, 0, 2, 100, levels=4)
        @test U2[4, 4] ≈ exact atol=1e-10
    end

    @testset "Diagonal improves accuracy: sin(x) on [0,π]" begin
        f = sin
        exact = 2.0
        U, _ = romberg(f, 0, π, 10, levels=4)
        @test U[4, 4] ≈ exact atol=1e-3
        U2, _ = romberg(f, 0, π, 100, levels=4)
        @test U2[4, 4] ≈ exact atol=1e-10
    end

    @testset "Diagonal improves accuracy: exp(x) on [0,1]" begin
        f = exp
        exact = ℯ - 1
        U, _ = romberg(f, 0, 1, 10, levels=4)
        @test U[4, 4] ≈ exact atol=1e-3
        U2, _ = romberg(f, 0, 1, 100, levels=4)
        @test U2[4, 4] ≈ exact atol=1e-10
    end

    # ── 3. Эффективный порядок точности ──────────────────────

    @testset "Diagonal strictly improves: sin(x) on [0,π]" begin
        # sin не является полиномом — уточнения реально нужны
        f = sin
        exact = 2.0
        U, _ = romberg(f, 0, π, 10, levels=4)
        for s in 2:4
            @test abs(U[s, s] - exact) < abs(U[s-1, s-1] - exact)
        end
        # Недостаточно данных для оценки на первых строках
        f = x -> x^2
        _, p_eff = romberg(f, 0, 1, 10, levels=4)
        for l in 1:4
            @test p_eff[1, l] == 0.0 || isnan(p_eff[1, l])
            @test p_eff[2, l] == 0.0 || isnan(p_eff[2, l])
        end
    end

    @testset "p_eff converges to 2 for trapezoid (x^2)" begin
        f = x -> x^2
        _, p_eff = romberg(f, 0, 1, 10, levels=5)
        # Первый столбец — порядок трапеции p=2
        @test p_eff[3, 1] ≈ 2.0 atol=0.1
        @test p_eff[4, 1] ≈ 2.0 atol=0.1
    end

    @testset "p_eff converges to 2 for trapezoid (sin)" begin
        f = sin
        _, p_eff = romberg(f, 0, π, 10, levels=5)
        @test p_eff[3, 1] ≈ 2.0 atol=0.2
        @test p_eff[4, 1] ≈ 2.0 atol=0.2
    end

    @testset "p_eff grows along columns" begin
        # На каждом следующем столбце порядок выше
        f = x -> x^4
        _, p_eff = romberg(f, 0, 1, 10, levels=5)
        @test p_eff[4, 2] > p_eff[4, 1]
    end

    # ── 4. Граничные случаи ───────────────────────────────────

    @testset "Constant function f(x) = 3" begin
        f = x -> 3.0
        exact = 9.0   # ∫₀³ 3 dx
        U, _ = romberg(f, 0, 3, 10, levels=4)
        # Трапеция точна для константы — весь первый столбец точный
        for s in 1:4
            @test U[s, 1] ≈ exact atol=1e-10
        end
        @test U[4, 4] ≈ exact atol=1e-10
    end

    @testset "Linear f(x) = 2x + 1 on [0, 2]" begin
        f = x -> 2x + 1
        exact = 6.0   # ∫₀² (2x+1) dx = [x²+x]₀² = 6
        U, _ = romberg(f, 0, 2, 10, levels=4)
        # Трапеция точна для линейных — U[s,1] точный при любом s
        for s in 1:4
            @test U[s, 1] ≈ exact atol=1e-10
        end
    end

    @testset "levels=1 returns 1x1 matrix" begin
        f = x -> x^2
        U, p_eff = romberg(f, 0, 1, 100, levels=1)
        @test size(U)     == (1, 1)
        @test size(p_eff) == (1, 1)
        @test U[1, 1] ≈ 1/3 atol=1e-4
    end

    @testset "r=3 (triple refinement)" begin
        f = x -> x^2
        exact = 1/3
        U, _ = romberg(f, 0, 1, 5, r=3, levels=4)
        @test U[4, 4] ≈ exact atol=1e-4
        U2, _ = romberg(f, 0, 1, 50, r=3, levels=4)
        @test U2[4, 4] ≈ exact atol=1e-10
    end

end