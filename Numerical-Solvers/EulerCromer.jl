function EulerCromer(x̄ₙ,p,Δt)
    # Euler Cromer for Lorenz system
    # 
    θ, μ, β = p
    xₙ, yₙ, zₙ = x̄ₙ

    ẋₙ = θ * (yₙ - xₙ)
    xₙ₊₁ = xₙ + ẋₙ * Δt
    
    ẏₙ = xₙ₊₁ * (μ - zₙ) - yₙ
    yₙ₊₁ =  yₙ + ẏₙ * Δt

    żₙ = xₙ₊₁ * yₙ₊₁ - β * zₙ
    zₙ₊₁ = zₙ + żₙ * Δt

    return [xₙ₊₁; yₙ₊₁; zₙ₊₁]
end
