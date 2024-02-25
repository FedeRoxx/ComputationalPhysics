using LinearAlgebra
# import Pkg; Pkg.add("Plots");
using Plots
using SpecialFunctions

function solution_FDM(N,dx)
    H = SymTridiagonal(2*ones(N), -ones(N-1)) / (dx^2)
    display(H)
    lambda_vec = eigvals(H)
    @show psi_vec = eigvecs(H)
    return lambda_vec, psi_vec
end

function solution_exact(N, dx)
    exact_lambda = [(pi*i)^2 for i in 1:N]
    psi_exact = zeros(N,N)
    for n in 1:N
        psi_exact[:,n] = [sqrt(2)*sin(n*i*pi*dx) for i in 1:N]
    end
    return exact_lambda, psi_exact
end

function L2_integral(v, dx)
    return sqrt(dx*sum([v_i^2 for v_i in v]))
end


N = 10
dx = 1/N

n_vec = [i for i in 1:N]
lambda_vec, psi_vec = solution_FDM(N,dx)
exact_lambda, psi_exact = solution_exact(N, dx)

# Plotting first 3 eigenfuntions
labels=["Numerical ψ1" "Numerical ψ2" "Numerical ψ3"]
plot(n_vec*dx, psi_vec[:,1:3], label=labels)
labels=["Exact ψ1" "Exact ψ2" "Exact ψ3"]
display(plot!(n_vec*dx, psi_exact[:,1:3], label=labels))

println(L2_integral(psi_exact[:,1], dx))
println(L2_integral(psi_vec[:,1],dx))



# plotting eigenval in function of n
plot(n_vec, exact_lambda, label="exact", show=true)
display(plot!(n_vec, lambda_vec, label="Numerical", show=true))



