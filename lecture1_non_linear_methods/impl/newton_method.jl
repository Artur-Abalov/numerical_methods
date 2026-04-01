using Printf, ForwardDiff

function balanced_coeff(ξ, τ = 1)
    return τ * ξ
end


function newton_method_simple(f, x0, ε)
    max_iters = 50
    x_nm = nothing
    x_n = x0
    x_n1 = nothing

    for i = 1:max_iters
        if i > 1
            x_nm = x_n
            x_n = x_n1
        end


        val = f(x_n)
        gradient = ForwardDiff.derivative(f, x_n)

        if abs(val) < 1e-14
            @printf "root found %.2f %d iters \n" x_n 0
            return x_n
        end

        if abs(gradient) < 1e-14 # попали в критическую точку, беда
            @printf "METHOD DIVERGED"
            return nothing
        end

        x_n1 = x_n - balanced_coeff(val/gradient)

        if i > 1
            stop_crit = abs((x_n1 - x_n)/(1 - (x_n1 - x_n)/(x_n - x_nm)))
            if(stop_crit < ε)
                root = x_n1
                @printf "root found %.2f %d iters \n" root i
                return root
            end
        end
    end

    @printf "METHOD DIVERGED"
    return nothing
end

function newton_method_augmented(f, x0, ε, k = 0.1)
    max_iters = 150
    x_nm = nothing
    x_n = x0
    x_n1 = nothing


    for i = 1:max_iters
        val = f(x_n)
        gradient = ForwardDiff.derivative(f, x_n)

        if abs(val) < 1e-14
            @printf "root found %.2f %d iters \n" x_n 0
            return x_n
        end

        if abs(gradient) < 1e-14 # попали в критическую точку, беда
            @printf "METHOD DIVERGED GRAD = 0"
            return nothing
        end

        x_n1_temp = x_n - val/gradient

        ψ_0 = f(x_n)^2
        ψ_1 = f(x_n1_temp)^2


        τ = (ψ_0 + k*ψ_1) / (ψ_0 + ψ_1)
        x_n1 = x_n - τ * val/gradient

        if x_nm != nothing
            stop_crit = abs((x_n1 - x_n)/(1 - (x_n1 - x_n)/(x_n - x_nm)))
            if(stop_crit < ε)
                root = x_n1
                @printf "root found %.2f %d iters \n" root i
                return root
            end
        end

        x_nm = x_n
        x_n = x_n1
    end

    @printf "METHOD DIVERGED"
    return nothing
end