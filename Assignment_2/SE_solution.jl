using LinearAlgebra
# import Pkg; Pkg.add("Plots");
using Plots
using SpecialFunctions

#####################
# PARTICLE IN A BOX #
#####################

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
    lambda_exact = [(pi*i)^2 for i in 1:N]
    psi_exact = zeros(N,N)
    for n in 1:N
        psi_exact[:,n] = [sqrt(2)*sin(n*i*pi*dx) for i in 1:N]
    end
    return lambda_exact, add_boundaries(psi_exact, N+1)
end

function L2_integral(v, dx)
    return sqrt(dx*sum([v_i*conj(v_i) for v_i in v]))
end

function add_boundaries(mat, N)
    # Adds zeros all around the matrix
    new_mat = zeros(N+1,N+1)
    new_mat[2:N,2:N] = mat
    return new_mat
end

# Initializing values

N = 100
dx = 1/N

x_vec = dx*[i for i in 0:N]
n_vec = [i for i in 1:N-1]
lambda_vec, psi_vec = solution_FDM(N-1,dx)
lambda_exact, psi_exact = solution_exact(N-1, dx)


# Plotting first 3 eigenfuntions
plotting=[2,3,4]
labels=["Numerical ψ1" "Numerical ψ2" "Numerical ψ3"]
plot(x_vec, psi_vec[:,plotting], label=labels)
labels=["Exact ψ1" "Exact ψ2" "Exact ψ3"]
display(plot!(x_vec, psi_exact[:,plotting], label=labels))
# display(plot(x_vec,  psi_exact[:,2].* psi_exact[:,4], label=labels))


println(L2_integral(psi_exact[:,1], dx))
println(L2_integral(psi_vec[:,1],dx))

#This to fix the sign
println(min(L2_integral(psi_exact[:,1]-psi_vec[:,1], dx), L2_integral(psi_exact[:,1]+psi_vec[:,1], dx)))
println(min(L2_integral(psi_exact[:,2]-psi_vec[:,2], dx), L2_integral(psi_exact[:,2]+psi_vec[:,2], dx)))
println(min(L2_integral(psi_exact[:,3]-psi_vec[:,3], dx), L2_integral(psi_exact[:,3]+psi_vec[:,3], dx)))

# Plotting eigenval in function of n
plot(n_vec, lambda_exact, label="Exact", show=true)
display(plot!(n_vec, lambda_vec, label="Numerical", show=true))




################################
# EXPANSIONS IN EIGENFUNCTIONS #
################################

# Defining overlaps integrals

function integrate_product(f1, f2, N, dx)
    # Returns integral of f1 * f2 
    return trapezoidal_integral(f1 .* f2, N, dx)
end

function trapezoidal_integral(f, N, dx)
    return (sum(f)-0.5*(f[1]-f[end])) * dx
end

function simpsons_integral(f, N, dx)
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

function initialize_from_psi1(N, dx)
    # This is eigenfunction ψ_1 
    psi_0 = [sqrt(2)*sin(i*pi*dx) for i in 0:N]
    return psi_0
end

function initialize_from_delta(N, dx)
    # This is a delta centered in 0.5
    if N % 2 != 0
        error("Number of intervals must be even.")
    end
    psi_0 = zeros(N+1)
    psi_0[Int(N/2)+1] = 1 / sqrt(dx)
    return psi_0
end

function time_evolution(alpha_n, psi_exact, lambda_exact, t, N)
    evoluted_Psi = zeros(N+1)im
    for n in 1:N-1
        evoluted_Psi[:] += alpha_n[n]*exp(-lambda_exact[n]*1im *t)*psi_exact[:,n+1]
    end
    return evoluted_Psi
end

# Testing overlaps between eigenfunctions
overlap = zeros(N-1,N-1)
for i in 2:N
    for j in 2:N
        overlap[i-1,j-1] = integrate_product(psi_exact[:,i], psi_exact[:,j], N, dx)
    end
end
println("<ψ_n|ψ_m>")
display(overlap)

# Time evolution from ψ_1
psi_0 = initialize_from_psi1(N, dx)
println("Norm of Ψ_0 as ψ_1: ", L2_integral(psi_0, dx))
alpha_n = [integrate_product(psi_exact[:,i], psi_0, N, dx) for i in 2:N]
println("α_n: ",  alpha_n)
evoluted_Psi = time_evolution(alpha_n, psi_exact, lambda_exact, 1, N)
println("Norm of Ψ(t)", L2_integral(evoluted_Psi, dx))

# Time evolution from δ(x' - 0.5)
initialize_from_delta(N, dx)
println("Norm of Ψ_0 as delta: ", L2_integral(psi_0, dx))
alpha_n = [integrate_product(psi_exact[:,i], psi_0, N, dx) for i in 2:N]
println("α_n: ",  alpha_n)
@show evoluted_Psi = time_evolution(alpha_n, psi_exact, lambda_exact, 1, N)
println("Norm of Ψ(t)", L2_integral(evoluted_Psi, dx))
# I have no idea what happens for t>0


##############################
# BOX WITH POTENTIAL BARRIER #
##############################

function solution_FDM_box(V_n, N, dx)
    H = SymTridiagonal(2*ones(N), -ones(N-1)) / (dx^2)
    H = H + diagm(V_n)
    display(H)
    lambda_vec = eigvals(H)
    @show psi_vec = eigvecs(H)
    for n in 1:N
        psi_vec[:,n] = psi_vec[:,n] / L2_integral(psi_vec[:,n], dx)
    end
    return lambda_vec, add_boundaries(psi_vec, N+1)
end

function potential_barrier(V0, N, dx)
    V_n = zeros(N-1)
    for i in 1:N-1
        if 1/3 <= i*dx < 2/3
            V_n[i] = V0
        end
    end
    return V_n
end

function initialize_from_superposition(psi_barrier)
    # This is eigenfunction ψ_1 
    psi_0 = (psi_barrier[:,2]+psi_barrier[:,3])/sqrt(2)
    return psi_0
end

function root_function(λ, v0)
    if λ > v0
        throw("Root function is defined for 0 < λ , v0")
    end
    k = sqrt(λ)
    κ = sqrt(v0 - λ)
    A = κ*sin(k/3)+k*cos(k/3)
    A = A^2 * exp(κ/3)
    B = κ*sin(k/3)-k*cos(k/3)
    B = B^2 * exp(-κ/3)
    return A - B
end

# Checking the results
V_n = zeros(N-1)
lambda_barrier, psi_barrier = solution_FDM_box(V_n, N-1, dx)

v0 = 22
V_n = potential_barrier(v0, N, dx)
lambda_barrier, psi_barrier = solution_FDM_box(V_n, N-1, dx)
# Plotting first 3 eigenfunctions
plotting=[1,2,3]
labels=["Numerical ψ1" "Numerical ψ2" "Numerical ψ3"]
display(plot(x_vec, psi_barrier[:,plotting], label=labels))
println(lambda_barrier)
#### NOTE ####
# The eigenvalues grow up in couples 
psi_0 = initialize_from_superposition(psi_barrier)
display(plot(x_vec, psi_0, title="Initial condition"))
println("Norm of Ψ_0 as (ψ_1 + ψ_2) / √2: ", L2_integral(psi_0, dx))
alpha_n = [integrate_product(psi_barrier[:,i], psi_0, N, dx) for i in 2:N]
println("α_n: ",  alpha_n)
evoluted_Psi = time_evolution(alpha_n, psi_barrier, lambda_barrier, 1, N)
println("Norm of Ψ(t)", L2_integral(evoluted_Psi, dx))

t_bar=π/(lambda_barrier[2]-lambda_barrier[1])
println(t_bar)

if false
anim = @animate for k in 1:200
    evoluted_Psi = time_evolution(alpha_n, psi_barrier, lambda_barrier, t_bar*k/100, N) #t_bar*(1-0.01*(k-50))
    module_Psi = [conj(psi_i)*psi_i for psi_i in evoluted_Psi]
    plot(x_vec, real.(module_Psi), ylim=(-0.05, 6), label="Squared module", xlabel="Index", ylabel="|Ψ|^2", title="Squared module at time τ*"*string(k/100))
end
display(gif(anim, "/home/frossi/ComputationalPhysics/Assignment_2/Time_evolution_barrier.gif", fps=60))
end
### It starts on one side and it end up on the other

#### Root finding
root_f = [root_function(λ_i, v0) for λ_i in 1:v0-1]
display(plot( 1:v0-1, root_f, xlabel="λ", label="f(λ)", title="Root function"))

# Here the root_finding algorithm, if I had one
#2/3 x cos(x/3) (6 sin(x/3) + x cos(x/3)) = 0
