function RK4(ẋ::Function, p, xₙ, Δt)

    k₁ = ẋ(p, xₙ)
    k₂ = ẋ(p, xₙ + 0.5*Δt*k₁)
    k₃ = ẋ(p, xₙ + 0.5*Δt*k₂)
    k₄ = ẋ(p, xₙ + Δt*k₃)

    xₙ₊₁ = xₙ + 1/6*Δt*(k₁ + 2*k₂ + 2*k₃ + k₄)

    return xₙ₊₁

end