function LorentzSystem(p,x̄)
    # Lorentz system
    # Returns the derivative of the state vector x̄
    θ=p[1]; μ=p[2]; β=p[3]
    x=x̄[1]; y=x̄[2]; z=x̄[3]
    return [θ*(y-x); 
            x*(μ-z)-y;
            x*y-β*z]
end
