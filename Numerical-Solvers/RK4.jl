include("lorentz_system.jl")

function RK4(p, xₙ, Δt,ẋ=LorentzSystem::Function)

    k₁ = ẋ(xₙ,p)
    k₂ = ẋ(xₙ + 0.5*Δt*k₁,p)
    k₃ = ẋ(xₙ + 0.5*Δt*k₂,p)
    k₄ = ẋ(xₙ + Δt*k₃,p)

    xₙ₊₁ = xₙ + 1/6*Δt*(k₁ + 2*k₂ + 2*k₃ + k₄)

    return xₙ₊₁

end