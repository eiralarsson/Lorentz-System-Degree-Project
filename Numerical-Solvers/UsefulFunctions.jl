function SolutionLorenz(p, Δt, x₀, N, Solver)
    X = zeros(3,N+1)
    X[:,1] = x₀
    x = x₀
    for i = 2:N+1
        X[:,i] = Solver(p, X[:,i-1], Δt)
    end
    return X
end

# Center has to be a column vector
function Points_on_Sphere(nr_of_dots, r, center)
    dots = zeros(3, nr_of_dots)
    for i=1:nr_of_dots
        dot = rand(-1.0:0.001:1.0,(1,3))
        dot = r*dot./norm(dot)
        dots[:,i] = dot + center
    end
    return dots
end

# Returns the solutions to lorenz with parameters p for all the initial 
# conditions in 'points' for one of our own solvers.
function PointsSolutions(p,Δt,N,Solver,points)
    solutions = []
    M = size(points)[2]
    for i = 1:M
        x = points[:,i] + displace
        X = SolutionLorenz(p, Δt, x, N, Solver)
        push!(solutions, X)
    end
    return solutions
end

# Same as above but with solvers from 'DifferentialEquations'
function PointsSolutionsJuliaSolver(p,Δt,N,Solver,points)
    solutions = []
    M = size(points)[2]
    for i = 1:M
        x = points[:,i]
        prob = ODEProblem(LorentzSystem,x,(0,Δt*N),[θ,μ,β])
        X = solve(prob, Solver())
        push!(solutions, X)
    end
    return solutions
end