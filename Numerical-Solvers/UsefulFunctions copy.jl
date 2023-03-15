using LinearAlgebra

function LorenzSolutionFixedTimeStep(p, Δt, N, x₀, Solver)
    if typeof(Solver) == UnionAll
        prob = ODEProblem(LorentzSystem,x₀,(0,Δt*N),p);
        X = solve(prob, Solver(), adaptive=false, dt=Δt);
        return X
    else
        X = zeros(3,N+1)
        X[:,1] = x₀
        x = x₀
        for i = 2:N+1
            X[:,i] = Solver(p, X[:,i-1], Δt)
        end
        return X
    end
end

function Correlation(x̄₁,x̄₂)
    # Returns the correlation between two vectors
    if x̄₁ == [0,0,0] && x̄₂== [0,0,0]
        return 1.0
    elseif x̄₁ == [0 0 0] && x̄₂== [0 0 0]
        return 1.0
    else
        return sum(x̄₁.*x̄₂)/norm(x̄₁)/norm(x̄₂)
    end
end

function CorrelationMatrix(M₁, M₂, num_init_points, N)
    M = zeros(num_init_points,N+1)
    for i = 1:num_init_points
        for j = 1:N+1
            M[i,j] = Correlation(M₁[3*(i-1)+1:3*i,j], M₂[3*(i-1)+1:3*i,j])
        end
    end
    return M
end

function PointsSolutions(p,Δt,N::Int64,initial_points,Solver, steps_per_Δt::Int64)
    solutions = zeros(3*length(initial_points), N+1)
    for i = eachindex(initial_points)
        x₀ = initial_points[i]
        X = LorenzSolutionFixedTimeStep(p, Δt/steps_per_Δt, N*steps_per_Δt, x₀, Solver)[:,1:steps_per_Δt:end]
        solutions[3*(i-1)+1:3*i,:] = X
    end
    return solutions
end

function EnergyFunction(x)
    # a very simple energy function
    return sum(x.^2) 
end