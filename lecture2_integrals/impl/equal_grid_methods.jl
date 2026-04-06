function create_grid(a, b, N)
   return a:(b-a)/N:b
end


function left_rectangle_method(f, a, b, N)
    if a == b return 0 end

    grid_points = create_grid(a, b, N)
    integral = 0
    for index in 1:length(grid_points)-1
        left = grid_points[index]

        v = f(left)
        integral += v * (b-a)/N
    end
    return integral

end

function right_rectangle_method(f, a, b, N)
    if a == b return 0 end
    grid_points = create_grid(a, b, N)
    integral = 0
    for index in 1:length(grid_points)-1
        right = grid_points[index+1]

        v = f(right)
        integral += v * (b-a)/N
    end
    return integral

end

function trapecion_method(f, a, b, N)
    if a == b return 0 end
    grid_points = create_grid(a, b, N)
    integral = 0

    for index in 1:length(grid_points)-1
        left = grid_points[index]
        right = grid_points[index+1]

        v = (f(left) + f(right))/2

        integral += v * (b-a)/N
    end
    return integral

end

function middle_point_method(f, a, b, N)
    if a == b return 0 end
    grid_points = create_grid(a, b, N)
    integral = 0

    for index in 1:length(grid_points)-1
        left = grid_points[index]
        right = grid_points[index+1]

        center = left + (right - left)/2

        v = f(center)

        integral += v * (b-a)/N
    end

    return integral
end