# numerical_methods

Учебный репозиторий с реализациями численных методов на Julia.

## Структура

- `lecture1_non_linear_methods`
  - `impl/`: методы решения нелинейных уравнений
    - дихотомия
    - метод Ньютона (обычный и модифицированный)
    - метод простых итераций
  - `test/`: тестовые и демонстрационные сценарии
- `lecture2_integrals`
  - `impl/`: численное интегрирование на равномерной сетке
    - левые/правые прямоугольники
    - трапеции
    - средняя точка
  - `test/`: тесты точности и сходимости
- `lecture3_runge_romberg`
  - `impl/`: уточнение интегралов Рунге–Ромберга
  - `test/`: тесты структуры, точности и порядка

## Требования

- Julia (1.x)
- Пакет `ForwardDiff` (для тестов метода Ньютона)

Установка зависимости:

```julia
import Pkg
Pkg.add("ForwardDiff")
```

## Запуск тестов

Из корня репозитория:

```bash
julia lecture1_non_linear_methods/test/newton_method.jl
julia lecture1_non_linear_methods/test/dihotomy_method.jl
julia lecture1_non_linear_methods/test/fixed_point_method.jl
julia lecture2_integrals/test/equal_grid_methods.jl
julia lecture3_runge_romberg/test/runge_romberg.jl
```
