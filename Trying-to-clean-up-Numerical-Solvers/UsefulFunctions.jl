using LinearAlgebra

function lorenz_solution_fixed_timestep(p, Δt, N, x₀, Solver)
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

function correlation(x̄₁,x̄₂)
    # Returns the correlation between two vectors
    if x̄₁ == [0,0,0] && x̄₂== [0,0,0]
        return 1.0
    else
        return sum(x̄₁.*x̄₂)/norm(x̄₁)/norm(x̄₂)
    end
end

function correlation_matrix(M₁, M₂)
    r = Integer(size(M₁)[1]/3)
    c = size(M₁)[2]
    M = zeros(r,c)
    for i = 1:r
        for j = 1:c
            M[i,j] = correlation(M₁[3*(i-1)+1:3*i,j], M₂[3*(i-1)+1:3*i,j])
        end
    end
    return M
end

function points_solutions_matrix(p,Δt,N::Int64,initial_points,solver; steps_per_Δt=1::Int64)
    solutions = zeros(3*size(initial_points)[2], N+1)
    for i = 1:size(initial_points)[2]
        x₀ = initial_points[:,i]
        X = lorenz_solution_fixed_timestep(p, Δt/steps_per_Δt, N*steps_per_Δt, x₀, solver)[:,1:steps_per_Δt:end]
        solutions[3*(i-1)+1:3*i,:] = X
    end
    return solutions
end

function energy_function(x)
    # a very simple energy function
    return sum(x.^2) 
end