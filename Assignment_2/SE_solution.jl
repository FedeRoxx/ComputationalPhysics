using LinearAlgebra
# import Pkg; Pkg.add("Plots");
using Plots
using SpecialFunctions

function solution_FDM(N,dx)
    H = SymTridiagonal(2*ones(N), -ones(N-1)) / (dx^2)
    display(H)
    lambda_vec = eigvals(H)
    @show psi_vec = eigvecs(H)
    for n in 1:N
        psi_vec[:,n] = psi_vec[:,n] / L2_integral(psi_vec[:,n], dx)
    end
    return lambda_vec, add_boundaries(psi_vec, N+1)
end

function solution_exact(N, dx)
    exact_lambda = [(pi*i)^2 for i in 1:N]
    psi_exact = zeros(N,N)
    for n in 1:N
        psi_exact[:,n] = [sqrt(2)*sin(n*i*pi*dx) for i in 1:N]
    end
    return exact_lambda, add_boundaries(psi_exact, N+1)
end

function L2_integral(v, dx)
    return sqrt(dx*sum([v_i^2 for v_i in v]))
end

function add_boundaries(mat, N)
    new_mat = zeros(N+1,N+1)
    new_mat[2:N,2:N] = mat
    return new_mat
end

function integrate_product(f1, f2, N, dx)
    return trapezoidal_integral(f1 .* f2, N, dx)
end

function trapezoidal_integral(f, N, dx)
    return (sum(f)-0.5*(f[1]-f[end])) * dx
end

function simpsons_method(f, N, dx)
    if N % 2 != 0
        error("Number of intervals must be even.")
    end

    sum = f[1]+f[N+1]
    
    for i in 2:N
        if i % 2 == 0
            sum += 2 * f[i]
        else
            sum += 4 * f[i]
        end
    end
    
    return sum * dx / 3
end



N = 10
dx = 1/N

x_vec = dx*[i for i in 0:N]
n_vec = [i for i in 1:N-1]
lambda_vec, psi_vec = solution_FDM(N-1,dx)
exact_lambda, psi_exact = solution_exact(N-1, dx)


# Plotting first 3 eigenfuntions
plotting=[2,3,4]
labels=["Numerical ψ1" "Numerical ψ2" "Numerical ψ3"]
plot(x_vec, psi_vec[:,plotting], label=labels)
labels=["Exact ψ1" "Exact ψ2" "Exact ψ3"]
display(plot!(x_vec, psi_exact[:,plotting], label=labels))
# display(plot(x_vec,  psi_exact[:,2].* psi_exact[:,4], label=labels))


println(L2_integral(psi_exact[:,1], dx))
println(L2_integral(psi_vec[:,1],dx))
println(min(L2_integral(psi_exact[:,1]-psi_vec[:,1], dx), L2_integral(psi_exact[:,1]+psi_vec[:,1], dx)))
println(min(L2_integral(psi_exact[:,2]-psi_vec[:,2], dx), L2_integral(psi_exact[:,2]+psi_vec[:,2], dx)))
println(min(L2_integral(psi_exact[:,3]-psi_vec[:,3], dx), L2_integral(psi_exact[:,3]+psi_vec[:,3], dx)))

# plotting eigenval in function of n
plot(n_vec, exact_lambda, label="exact", show=true)
display(plot!(n_vec, lambda_vec, label="Numerical", show=true))

overlap = zeros(N-1,N-1)
for i in 2:N
    for j in 2:N
        overlap[i-1,j-1] = integrate_product(psi_exact[:,i], psi_exact[:,j], N, dx)
    end
end
display(overlap)
println(psi_exact[:,6], psi_exact[:,7])
integrate_product(psi_exact[:,6], psi_exact[:,7], N, dx)





