using LinearAlgebra

function Beta(x, θ, μ)
    return (x[1]*x[2]*(θ+μ)-θ*(x[1]^2)-x[2]^2)/(x[3]^2)
end

function Mu(x, θ, β)
    return θ*x[1]/x[2] + x[2]/x[1] + β*x[3]^2/x[2]/x[1]
end

jahY71&#218kAA8 ö