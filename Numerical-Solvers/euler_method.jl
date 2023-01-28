include("lorentz_system.jl")


function EulerForward(p,x̄,Δt,ẋ=LorentzSystem::Function)
    x̄ = x̄ + Δt*ẋ(x̄,p)
    return x̄
end

