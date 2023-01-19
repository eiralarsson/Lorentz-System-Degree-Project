function LorentzSystem(sigma,rho,beta,x,y,z)
    return [sigma*(y-x); 
            x*(rho-z)-y;
            x*y-beta*z]
end
