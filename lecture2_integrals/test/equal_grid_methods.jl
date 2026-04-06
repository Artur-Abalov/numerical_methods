include("../impl/equal_grid_methods.jl")

using Test

function all_methods(f, a, b, N)
    return (
        left  = left_rectangle_method(f, a, b, N),
        right = right_rectangle_method(f, a, b, N),
        trap  = trapecion_method(f, a, b, N),
        mid   = middle_point_method(f, a, b, N),
    )
end

# ─────────────────────────────────────────────────────────────
#  Тесты
# ─────────────────────────────────────────────────────────────

@testset "Numerical Integration" begin

    # ── 1. Базовые случаи ─────────────────────────────────────

    @testset "Constant function f(x) = 5" begin
        f = x -> 5.0
        r = all_methods(f, 0, 4, 100)
        # Все методы точны для константы
        @test r.left  ≈ 20.0 atol=1e-10
        @test r.right ≈ 20.0 atol=1e-10
        @test r.trap  ≈ 20.0 atol=1e-10
        @test r.mid   ≈ 20.0 atol=1e-10
    end

    @testset "Linear f(x) = x on [0, 1]" begin
        f = x -> x
        exact = 0.5
        r = all_methods(f, 0, 1, 1000)
        # Трапеция и средняя точка точны для линейных функций
        @test r.trap ≈ exact atol=1e-10
        @test r.mid  ≈ exact atol=1e-10
        # Прямоугольники сходятся, но с погрешностью O(h)
        @test abs(r.left  - exact) < 1e-3
        @test abs(r.right - exact) < 1e-3
    end

    @testset "Linear f(x) = 2x + 3 on [1, 5]" begin
        f = x -> 2x + 3
        exact = 36.0   # ∫₁⁵ (2x+3) dx = [x²+3x]₁⁵ = 40−4 = 36
        r = all_methods(f, 1, 5, 500)
        @test r.trap ≈ exact atol=1e-10
        @test r.mid  ≈ exact atol=1e-10
        # Для левого/правого прямоугольника: err ≈ h*(f(b)-f(a))/2
        # h = 4/500 = 0.008, f(b)-f(a) = 8 → err ≈ 0.032
        @test r.left  ≈ exact atol=0.05
        @test r.right ≈ exact atol=0.05
    end

    # ── 2. Полиномы ───────────────────────────────────────────

    @testset "Quadratic f(x) = x^2 on [0, 1]" begin
        f = x -> x^2
        exact = 1/3
        r = all_methods(f, 0, 1, 1000)
        @test r.left  ≈ exact atol=1e-3
        @test r.right ≈ exact atol=1e-3
        @test r.trap  ≈ exact atol=1e-5
        @test r.mid   ≈ exact atol=1e-5
    end

    @testset "Cubic f(x) = x^3 on [0, 2]" begin
        f = x -> x^3
        exact = 4.0    # ∫₀² x³ dx = [x⁴/4]₀² = 4
        r = all_methods(f, 0, 2, 1000)
        @test r.left  ≈ exact atol=1e-2
        @test r.right ≈ exact atol=1e-2
        @test r.trap  ≈ exact atol=1e-5
        @test r.mid   ≈ exact atol=1e-5
    end

    @testset "Polynomial f(x) = x^4 - 2x + 1 on [-1, 1]" begin
        f = x -> x^4 - 2x + 1
        exact = 2.4    # ∫₋₁¹ (x⁴−2x+1) dx = [x⁵/5 − x² + x]₋₁¹ = 12/5
        r = all_methods(f, -1, 1, 1000)
        @test r.left  ≈ exact atol=1e-2
        @test r.right ≈ exact atol=1e-2
        @test r.trap  ≈ exact atol=1e-4
        @test r.mid   ≈ exact atol=1e-4
    end

    # ── 3. Классические функции ───────────────────────────────

    @testset "sin(x) on [0, π]" begin
        f = sin
        exact = 2.0
        r = all_methods(f, 0, π, 1000)
        @test r.left  ≈ exact atol=1e-4
        @test r.right ≈ exact atol=1e-4
        @test r.trap  ≈ exact atol=1e-5
        @test r.mid   ≈ exact atol=1e-5
        # На [0,π] sin симметричен — погрешности left и right близки
        @test abs(r.left - exact) ≈ abs(r.right - exact) atol=1e-5
    end

    @testset "cos(x) on [0, π/2]" begin
        f = cos
        exact = 1.0
        r = all_methods(f, 0, π/2, 1000)
        @test r.left  ≈ exact atol=1e-3
        @test r.right ≈ exact atol=1e-3
        @test r.trap  ≈ exact atol=1e-5
        @test r.mid   ≈ exact atol=1e-5
    end

    @testset "cos(x) on [0, 2π] = 0" begin
        f = cos
        r = all_methods(f, 0, 2π, 1000)
        @test r.left  ≈ 0.0 atol=1e-13
        @test r.right ≈ 0.0 atol=1e-13
        @test r.trap  ≈ 0.0 atol=1e-13
        @test r.mid   ≈ 0.0 atol=1e-13
    end

    @testset "exp(x) on [0, 1]" begin
        f = exp
        exact = ℯ - 1        # ≈ 1.71828…
        r = all_methods(f, 0, 1, 1000)
        @test r.left  ≈ exact atol=1e-3
        @test r.right ≈ exact atol=1e-3
        @test r.trap  ≈ exact atol=1e-5
        @test r.mid   ≈ exact atol=1e-5
    end

    @testset "1/x on [1, e]" begin
        f = x -> 1/x
        exact = 1.0       # ∫₁ᵉ 1/x dx = ln(e) = 1
        r = all_methods(f, 1, ℯ, 1000)
        @test r.left  ≈ exact atol=1e-3
        @test r.right ≈ exact atol=1e-3
        @test r.trap  ≈ exact atol=1e-5
        @test r.mid   ≈ exact atol=1e-5
    end

    # ── 4. Сходимость — порядок точности ─────────────────────

    @testset "Convergence order: left/right = O(h)" begin
        f = x -> x^2
        a, b = 0.0, 1.0
        exact = 1/3
        err10   = abs(left_rectangle_method(f, a, b, 10)   - exact)
        err100  = abs(left_rectangle_method(f, a, b, 100)  - exact)
        err1000 = abs(left_rectangle_method(f, a, b, 1000) - exact)
        # При N×10 ошибка должна падать в ~10 раз
        @test err10  / err100  ≈ 10.0 atol=1.5
        @test err100 / err1000 ≈ 10.0 atol=1.0
    end

    @testset "Convergence order: trapezoid = O(h^2)" begin
        f = x -> x^2
        a, b = 0.0, 1.0
        exact = 1/3
        err10   = abs(trapecion_method(f, a, b, 10)   - exact)
        err100  = abs(trapecion_method(f, a, b, 100)  - exact)
        err1000 = abs(trapecion_method(f, a, b, 1000) - exact)
        # При N×10 ошибка падает в ~100 раз
        @test err10  / err100  ≈ 100.0 atol=20.0
        @test err100 / err1000 ≈ 100.0 atol=20.0
    end

    @testset "Convergence order: midpoint = O(h^2)" begin
        f = x -> x^2
        a, b = 0.0, 1.0
        exact = 1/3
        err10   = abs(middle_point_method(f, a, b, 10)   - exact)
        err100  = abs(middle_point_method(f, a, b, 100)  - exact)
        err1000 = abs(middle_point_method(f, a, b, 1000) - exact)
        @test err10  / err100  ≈ 100.0 atol=20.0
        @test err100 / err1000 ≈ 100.0 atol=20.0
    end

    @testset "Midpoint is more accurate than trapezoid (same N)" begin
        # Для гладких функций: ошибка средней точки ≈ 0.5 × ошибки трапеции
        f = x -> x^2
        a, b = 0.0, 1.0
        exact = 1/3
        err_trap = abs(trapecion_method(f, a, b, 10)    - exact)
        err_mid  = abs(middle_point_method(f, a, b, 10) - exact)
        @test err_mid < err_trap
    end

    # ── 5. Граничные и специальные случаи ─────────────────────

    @testset "N = 1 (minimum partition)" begin
        f = x -> x
        # Не должно бросать исключений; результат — грубая оценка
        @test_nowarn left_rectangle_method(f, 0, 1, 1)
        @test_nowarn right_rectangle_method(f, 0, 1, 1)
        @test_nowarn trapecion_method(f, 0, 1, 1)
        @test_nowarn middle_point_method(f, 0, 1, 1)
    end

    @testset "a == b (degenerate interval)" begin
        f = x -> x^2
        @test left_rectangle_method(f, 2, 2, 100)  ≈ 0.0 atol=1e-15
        @test right_rectangle_method(f, 2, 2, 100) ≈ 0.0 atol=1e-15
        @test trapecion_method(f, 2, 2, 100)        ≈ 0.0 atol=1e-15
        @test middle_point_method(f, 2, 2, 100)     ≈ 0.0 atol=1e-15
    end

    @testset "Reversed limits a > b" begin
        f = x -> x
        # ∫₁⁰ x dx = −∫₀¹ x dx = −0.5
        exact = -0.5
        @test left_rectangle_method(f, 1, 0, 1000)  ≈ exact atol=1e-3
        @test right_rectangle_method(f, 1, 0, 1000) ≈ exact atol=1e-3
        @test trapecion_method(f, 1, 0, 1000)        ≈ exact atol=1e-10
        @test middle_point_method(f, 1, 0, 1000)     ≈ exact atol=1e-10
    end

    @testset "Negative integrand f(x) = -x^2" begin
        f = x -> -x^2
        exact = -1/3
        r = all_methods(f, 0, 1, 1000)
        @test r.left  ≈ exact atol=1e-3
        @test r.right ≈ exact atol=1e-3
        @test r.trap  ≈ exact atol=1e-5
        @test r.mid   ≈ exact atol=1e-5
    end

    @testset "Sign-changing integrand f(x) = sin(x) on [0, 2π] = 0" begin
        f = sin
        r = all_methods(f, 0, 2π, 1000)
        @test r.left  ≈ 0.0 atol=1e-13
        @test r.right ≈ 0.0 atol=1e-13
        @test r.trap  ≈ 0.0 atol=1e-13
        @test r.mid   ≈ 0.0 atol=1e-13
    end

    @testset "Large interval f(x) = 1 on [0, 10000]" begin
        f = x -> 1.0
        exact = 10_000.0
        r = all_methods(f, 0, 10_000, 100)
        @test r.left  ≈ exact atol=1e-6
        @test r.right ≈ exact atol=1e-6
        @test r.trap  ≈ exact atol=1e-6
        @test r.mid   ≈ exact atol=1e-6
    end

    @testset "Negative interval f(x) = x on [-3, -1]" begin
        f = x -> x
        exact = -4.0   # ∫₋₃⁻¹ x dx = [x²/2]₋₃⁻¹ = 0.5 − 4.5 = −4
        r = all_methods(f, -3, -1, 1000)
        @test r.trap ≈ exact atol=1e-10
        @test r.mid  ≈ exact atol=1e-10
        @test r.left  ≈ exact atol=5e-3
        @test r.right ≈ exact atol=5e-3
    end

end

