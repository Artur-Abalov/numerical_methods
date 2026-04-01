using Printf

function fixed_point_method(φ, x0, ε)
    max_iters = 10000
    x_nm = nothing
    x_n = x0
    x_n1 = φ(x_n)
    for i = 1:max_iters


        x_nm = x_n
        x_n = x_n1
        x_n1 = φ(x_n)



        stop_crit = abs((x_n1 - x_n)/(1 - (x_n1 - x_n)/(x_n - x_nm)))

        if(stop_crit < ε)
            root = x_n
            @printf "root found %.2f %d iters \n" root i
            return root
        end

    end

    @printf "METHOD DIVERGED"
    return nothing
end

φ(x) = (x + 2)^(1/3)

fixed_point_method(φ, 1.5, 1e-5)