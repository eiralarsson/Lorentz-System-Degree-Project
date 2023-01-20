function EulerForward(ẋ::Function,p,x̄,Δt)
    x̄ = x̄ + Δt*ẋ(p,x̄)
    return x̄
end
