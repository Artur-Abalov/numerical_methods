using Printf

function dihotomy(f, x0, x1, ε)
    max_iters = 50

    if f(x0) == 0
        @printf "root found %.2f iter count %d \n" x0 0

       return x0
    end

    if f(x1) == 0
        @printf "root found %.2f iter count %d \n" x1 0

       return x1
    end


    if f(x0) * f(x1) > 0
        @printf "same signs"
        return nothing
    end

    x_n = x0
    x_n1 = x1

    for i = 1:max_iters

        x_c = (x_n + x_n1)/2

        if f(x_c) * f(x_n) < 0
            x_n1 = x_c
        else
            x_n = x_c
        end

        if abs(x_n1 - x_n) < ε
            root = (x_n + x_n1) / 2
            @printf "root found %.2f iter count %d \n" root i
            return root
        end
    end

end

function f(x)
    (x-1) * (x-5)
end

dihotomy(f, -3, 4, 0.1)