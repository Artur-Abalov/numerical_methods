include("../impl/newton_method.jl")

using Test

@testset "Newton Method Tests" begin

    @testset "Простые корни" begin
        # x² - 4 = 0, корни: x = ±2
        f1(x) = x^2 - 4
        @test newton_method_simple(f1, 1.0, 1e-10) ≈ 2.0  atol=1e-8
        @test newton_method_simple(f1, -1.0, 1e-10) ≈ -2.0 atol=1e-8

        # x³ - 8 = 0, корень: x = 2
        f2(x) = x^3 - 8
        @test newton_method_simple(f2, 1.0, 1e-10) ≈ 2.0  atol=1e-8
    end

    @testset "Тригонометрия" begin
        # sin(x) = 0, ближайший корень к 3.0 — это π
        f3(x) = sin(x)
        @test newton_method_simple(f3, 3.0, 1e-10) ≈ π  atol=1e-8

        # cos(x) = 0, ближайший корень к 1.0 — это π/2
        f4(x) = cos(x)
        @test newton_method_simple(f4, 1.0, 1e-10) ≈ π/2  atol=1e-8
    end

    @testset "Экспонента и логарифм" begin
        # exp(x) - 2 = 0, корень: x = ln(2)
        f5(x) = exp(x) - 2
        @test newton_method_simple(f5, 1.0, 1e-10) ≈ log(2)  atol=1e-8

        # log(x) - 1 = 0, корень: x = e
        f6(x) = log(x) - 1
        @test newton_method_simple(f6, 2.0, 1e-10) ≈ ℯ  atol=1e-8
    end

    @testset "Точность сходимости" begin
        f7(x) = x^2 - 2  # корень: √2
        
        root_low  = newton_method_simple(f7, 1.0, 1e-4)
        root_high = newton_method_simple(f7, 1.0, 1e-10)

        # Более жёсткий допуск — более точный результат
        @test abs(root_low  - sqrt(2)) < 1e-4
        @test abs(root_high - sqrt(2)) < 1e-10
    end

    @testset "Крайние случаи" begin
        # Начальное приближение уже является корнем
        f8(x) = x^2 - 9
        @test newton_method_simple(f8, 3.0, 1e-10) ≈ 3.0  atol=1e-8

        # Нет вещественных корней — метод должен вернуть nothing
        f9(x) = x^2 + 1
        @test newton_method_simple(f9, 1.0, 1e-10) === nothing
    end

end