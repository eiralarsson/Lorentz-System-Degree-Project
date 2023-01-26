function EulerForward(ẋ::Function,p,x̄,Δt)
    x̄ = x̄ + Δt*ẋ(x̄,p)
    return x̄
end

