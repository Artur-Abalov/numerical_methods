include("../impl/fixed_point_method.jl")

println("==== TEST 1: cubic ====")
φ1(x) = (x + 2)^(1/3)
fixed_point_method(φ1, 1.5, 1e-5)

println("\n==== TEST 2: cos(x) ====")
φ2(x) = cos(x)
fixed_point_method(φ2, 0.5, 1e-6)

println("\n==== TEST 3: linear ====")
φ3(x) = 0.5x + 1
fixed_point_method(φ3, 0.0, 1e-6)

println("\n==== TEST 4: slow convergence - беда сошлось при 1000+ итераций ====")
φ4(x) = 0.99x + 0.01
fixed_point_method(φ4, 0.0, 1e-6)

println("\n==== TEST 5: divergence ====")
φ5(x) = 2x + 1
fixed_point_method(φ5, 0.0, 1e-6)

println("\n==== TEST 6: oscillation - метод зациклился ====")
φ6(x) = -x
fixed_point_method(φ6, 1.0, 1e-6)

println("\n==== TEST 7: sqrt domain ====")
φ7(x) = sqrt(x + 2)
fixed_point_method(φ7, 1.0, 1e-6)
