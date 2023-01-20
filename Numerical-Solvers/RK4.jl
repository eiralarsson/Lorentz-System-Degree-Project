function RK4(ẋ::Function, xₙ, tₙ , Δt)

    k₁ = ẋ(tₙ, xₙ)
    k₂ = ẋ(tₙ + 0.5*Δt, xₙ + 0.5*Δt*k₁)
    k₃ = ẋ(tₙ + 0.5*Δt, xₙ + 0.5*Δt*k₂)
    k₄ = ẋ(tₙ + Δt, xₙ + Δt*k₃)

    xₙ₊₁ = xₙ + 1/6*Δt*(k₁ + 2*k₂ + 2*k₃ + k₄)
    

end