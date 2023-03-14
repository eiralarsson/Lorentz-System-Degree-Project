function LorentzSystem(x̄,P,t=1)
    # Lorentz system
    # Returns the derivative of the state vector x̄
    p = P(x̄)
    θ=p[1]; μ=p[2]; β=p[3]
    x=x̄[1]; y=x̄[2]; z=x̄[3]
    return [θ*(y-x); 
            x*(μ-z)-y;
            x*y-β*z]
end

function LorentzSystemBeta(x̄,P,t=1)
    # Lorentz system
    # Returns the derivative of the state vector x̄
    p = P(x̄)
    θ=p[1]; μ=p[2];
    x=x̄[1]; y=x̄[2]; z=x̄[3]
    a=p[3]; b=p[4]; c=p[5]
    β=(x*y - ((y - x*(μ - z))*(b - y) - θ*(a - x)*(x - y))/(c - z))/z
    return [θ*(y-x); 
            x*(μ-z)-y;
            x*y-β*z]
end