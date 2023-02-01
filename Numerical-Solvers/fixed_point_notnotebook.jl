using Plots
using DifferentialEquations
include("EulerCromer.jl")
fixed_point = [0,0,0]
x₀ = fixed_point
Δt = 0.001
N = 100
X̄_EulerCromer = zeros(3,N)
X̄_EulerCromer[:,1] = x₀ 
p = [10,8/3,23.5]

for i = 2:N
    X̄_EulerCromer[:,i] = EulerCromer(X̄_EulerCromer[:,i-1], p, Δt)
    print(i)
end
print(X̄_EulerCromer[:,N])