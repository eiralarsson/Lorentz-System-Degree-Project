using LinearAlgebra
using Plots
include("lorentz_system.jl")


function DTM(p,x̄,Δt, terms=10::Int64)
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
        # I honestly don't know why it works to use 2*Δt instead of Δt, but it does
        x̄[1] = x̄[1] + C[1,i]*(2*Δt)^(i-1)
        x̄[2] = x̄[2] + C[2,i]*(2*Δt)^(i-1)
        x̄[3] = x̄[3] + C[3,i]*(2*Δt)^(i-1)
    end

    return x̄
end


function EulerForward(p,x̄,Δt,ẋ=LorentzSystem::Function)
    x̄ = x̄ + Δt*ẋ(x̄,p)
    return x̄
end


function EulerCromer(p,x̄ₙ,Δt)
    # Euler Cromer for Lorenz system
    # 
    θ, μ, β = p
    xₙ, yₙ, zₙ = x̄ₙ

    ẋₙ = θ * (yₙ - xₙ)
    xₙ₊₁ = xₙ + ẋₙ * Δt
    
    ẏₙ = xₙ * (μ - zₙ) - yₙ
    yₙ₊₁ =  yₙ + ẏₙ * Δt

    żₙ = xₙ₊₁ * yₙ₊₁ - β * zₙ
    zₙ₊₁ = zₙ + żₙ * Δt

    return [xₙ₊₁; yₙ₊₁; zₙ₊₁]
end


function RK4(p, xₙ, Δt,ẋ=LorentzSystem::Function)

    k₁ = ẋ(xₙ,p)
    k₂ = ẋ(xₙ + 0.5*Δt*k₁,p)
    k₃ = ẋ(xₙ + 0.5*Δt*k₂,p)
    k₄ = ẋ(xₙ + Δt*k₃,p)

    xₙ₊₁ = xₙ + 1/6*Δt*(k₁ + 2*k₂ + 2*k₃ + k₄)

    return xₙ₊₁

end