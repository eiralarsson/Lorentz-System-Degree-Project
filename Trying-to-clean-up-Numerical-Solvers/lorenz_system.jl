function LorenzSystem(x̄,p,t=1)
    # Lorenz system
    # Returns the derivative of the state vector x̄
    θ, μ, β = p
    x=x̄[1]; y=x̄[2]; z=x̄[3]
    return [θ*(y-x); 
            x*(μ-z)-y;
            x*y-β*z]
end
