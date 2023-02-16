using LinearAlgebra

function SolutionLorenz(p, Δt, N::Int64, x₀, Solver)
    X = zeros(3,N+1)
    X[:,1] = x₀
    x = x₀
    for i = 2:N+1
        X[:,i] = Solver(p, X[:,i-1], Δt)
    end
    return X
end

function SolutionLorenzJuliaSolver(p, Δt, N::Int64, x₀, Solver)
    prob = ODEProblem(LorentzSystem,x₀,(0,Δt*N),p);
    X = solve(prob, Solver(), adaptive=false, dt=Δt);
    return X
end

# Center has to be a column vector
function Points_on_Sphere(nr_of_dots, r, center=[0,0,0]')
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
function PointsSolutions(p,Δt,N,initial_points,Solver)
    solutions = []
    M = size(initial_points)[2]
    for i = 1:M
        x = initial_points[:,i]
        X = SolutionLorenz(p, Δt, x, N, Solver)
        push!(solutions, X)
    end
    return solutions
end

# Same as above but with solvers from 'DifferentialEquations'
function PointsSolutionsJuliaSolver(p,Δt,N,initial_points,Solver)
    solutions = []
    M = size(initial_points)[2]
    for i = 1:M
        x = initial_points[:,i];
        prob = ODEProblem(LorentzSystem,x,(0,Δt*N),p);
        X = solve(prob, Solver(), adaptive=false, dt=Δt);
        push!(solutions, X)
    end
    return solutions
end

function Correlation(vec1,vec2)
    return sum(vec1.*vec2)/norm(vec1)/norm(vec2)
end

function CorrelationMatrix(p,N,dt,initial_points,Solvers,ratio)
    r = ratio
    M = zeros(length(initial_points), N+1)
    for i = eachindex(initial_points)
        x = initial_points[i]
        X = []
        for k = 1:2
            if typeof(Solvers[k]) == UnionAll
                push!(X,SolutionLorenzJuliaSolver(p, dt/r[k], N*r[k], x, Solvers[k])[1:ratio[k]:end])
            else
                push!(X,SolutionLorenz(p, dt/r[k], N*r[k], x, Solvers[k])[:,1:ratio[k]:end])
            end
        end
        for j = 1:N+1
            M[i,j] = Correlation(X[1][:,j], X[2][:,j])
        end
    end
    return M
end

function PointsInCuboid(dx,xint,dy,yint,dz,zint)
    vec = []
    for x = xint[1]:dx:xint[2]
        for y = yint[1]:dy:yint[2]
            for z = zint[1]:dz:zint[2]
                push!(vec,[x,y,z])
            end
        end
    end
    return vec
end