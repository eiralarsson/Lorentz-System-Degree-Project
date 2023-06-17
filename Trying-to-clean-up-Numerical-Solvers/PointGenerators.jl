function points_in_cuboid(I₁,I₂,I₃)
    vec = zeros(3,size(I₁)[1]*size(I₂)[1]*size(I₃)[1])
    i = 1
    for x₁ = I₁
        for x₂ = I₂
            for x₃ = I₃
                vec[:,i] = [x₁, x₂, x₃]
                i = i + 1
            end
        end
    end
    return vec
end

function points_on_line(x̄₀, v̄, Iₛ , Δs)
    vec = []
    for s = Iₛ[1]:Δs:Iₛ[2]
        x = x̄₀ + s * v̄
        push!(vec,x)
    end
    return vec
end

# Center has to be a column vector
function points_on_sphere(nr_of_dots, r; center=[0 0 0])
    dots = zeros(3, nr_of_dots)
    for i=1:nr_of_dots
        dot = rand(-1.0:0.001:1.0,(1,3))
        dot = r*dot./norm(dot)
        dots[:,i] = dot + center
    end
    return dots
end