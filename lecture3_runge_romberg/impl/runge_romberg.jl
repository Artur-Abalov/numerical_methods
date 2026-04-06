include("../../lecture2_integrals/impl/equal_grid_methods.jl")

function romberg(f, a, b, N; r=2, levels=4)
    m = levels
    U = zeros(m, m)
    R = zeros(m, m)
    q = 2 # trapecion
    p = 2
    p_eff = zeros(m, m)

    for i in 1:m
        U[i, 1] = trapecion_method(f, a, b, N * r^(i-1))
    end


    for s in 2:m
        for l in 1:s-1
            R[s, l] = (U[s, l] - U[s-1, l]) / (r^(p + (l-1) * q) - 1)
            U[s, l + 1] = U[s, l] + R[s,l]
        end
    end

    for s in 3:m
       for l in 1:s-2
           p_eff[s, l] = log(abs(R[s-1, l]/R[s,l]))/log(r)
       end
    end

    return U, p_eff
end

