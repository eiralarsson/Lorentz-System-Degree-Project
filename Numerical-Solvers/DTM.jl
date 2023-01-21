using LinearAlgebra
using Plots

# Implementations of the differential transformation method (DTM) for lorenz system with 10 terms
# The DTM is a method for solving differential equations with a finite number of terms
function DTM(p, x̄, terms::Int64,Δt)
    θ=p[1]; μ=p[2]; β=p[3]
    C = zeros(3, terms)
    C[1,1] = x̄[1]; C[2,1] = x̄[2]; C[3,1] = x̄[3]
    for i=2:terms
        C[1,i] = θ*(C[2,i-1]-C[1,i-1])/i
        C[2,i] = (C[1,i-1]*μ-sum(C[1,l]*C[3,i-l] for l in 1:i-1) - C[2,i-1])/i
        C[3,i] = (sum(C[1,l]*C[2,i-l] for l in 1:i-1) - β*C[3,i-1])/i
    end

    x̄ = zeros(3)
    for i=1:terms
        x̄[1] = x̄[1] + C[1,i]*Δt^(i-1)
        x̄[2] = x̄[2] + C[2,i]*Δt^(i-1)
        x̄[3] = x̄[3] + C[3,i]*Δt^(i-1)
    end

    return x̄
end

