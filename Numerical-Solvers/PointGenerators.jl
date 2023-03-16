
function PointsInCuboid(Δx₁,I₁,Δx₂,I₂,Δx₃,I₃)
    vec = []
    for x₁ = I₁[1]:Δx₁:I₁[2]
        for x₂ = I₂[1]:Δx₂:I₂[2]
            for x₃ = I₃[1]:Δx₃:I₃[2]
                push!(vec,[x₁,x₂,x₃])
            end
        end
    end
    return vec
end

function PointsOnLine(x̄₀, v̄, Iₛ , Δs)
    vec = []
    for s = Iₛ[1]:Δs:Iₛ[2]
        x = x̄₀ + s * v̄
        push!(vec,x)
    end
    return vec
end

# Center has to be a column vector
function PointsOnSphere(nr_of_dots, r, center=[0 0 0])
    dots = zeros(3, nr_of_dots)
    for i=1:nr_of_dots
        dot = rand(-1.0:0.001:1.0,(1,3))
        dot = r*dot./norm(dot)
        dots[:,i] = dot + center
    end
    return dots
end