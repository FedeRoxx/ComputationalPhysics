using LinearAlgebra
using Plots
using SpecialFunctions
using LaTeXStrings
using Printf

#####################
# PARTICLE IN A BOX #
#####################

function solution_FDM(N,dx)
    H = SymTridiagonal(2*ones(N), -ones(N-1)) / (dx^2)
    # display(H)
    lambda_vec = eigvals(H)
    psi_vec = eigvecs(H)
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
plot(x_vec, psi_vec[:,plotting], label=labels, xlabel="x", ylabel="ψ(x)", title="First 3 eigenfunctions")
labels=["Exact ψ1" "Exact ψ2" "Exact ψ3"]
display(plot!(x_vec, psi_exact[:,plotting], label=labels))
labels=[L"ψ_1" L"ψ_2" L"ψ_3"]
display(plot(x_vec, psi_exact[:,plotting], label=labels, xlabel=L"x", ylabel=L"ψ(x)", title="First 3 eigenfunctions "*L"ψ_1 \, ψ_2 \, ψ_3"))


println("L2 norm of psi_exact_1 : ", L2_integral(psi_exact[:,2], dx))
println("L2 norm of psi_numerical_1 : ", L2_integral(psi_vec[:,2],dx))

#This to fix the sign
# println(min(L2_integral(psi_exact[:,1]-psi_vec[:,1], dx), L2_integral(psi_exact[:,1]+psi_vec[:,1], dx)))
println("L2 norm of psi_exact_1 - psi_numerical_1 : ", min(L2_integral(psi_exact[:,2]-psi_vec[:,2], dx), L2_integral(psi_exact[:,2]+psi_vec[:,2], dx)))
println("L2 norm of psi_exact_2 - psi_numerical_2 : ", min(L2_integral(psi_exact[:,3]-psi_vec[:,3], dx), L2_integral(psi_exact[:,3]+psi_vec[:,3], dx)))

# Plotting eigenval in function of n
plot(n_vec, lambda_exact, xlabel=L"n", ylabel=L"λ_n", label=L"λ_n"*" exact", title="Comparison of eigenvalues as function of "*L"n")
display(plot!(n_vec, lambda_vec, label=L"λ_n"*" numerical"))
display(plot(n_vec, (lambda_exact.-lambda_vec)./lambda_exact, label=L"λ_n"*" relative difference"))


# Plotting numerical vs exact_solution
display(plot(lambda_exact, lambda_vec, label="Numerical"))





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

function initialize_from_psi1_exact(N, dx)
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
psi_0 = initialize_from_psi1_exact(N, dx)
println("Norm of Ψ_0 as ψ_1: ", L2_integral(psi_0, dx))
alpha_n = [integrate_product(psi_vec[:,i], psi_0, N, dx) for i in 2:N]
println("α_n: ",  alpha_n)
evoluted_Psi = time_evolution(alpha_n, psi_exact, lambda_exact, pi/lambda_exact[1], N)
println("Norm of Ψ(t): ", L2_integral(evoluted_Psi, dx))
module_Psi = [conj(psi_i)*psi_i for psi_i in evoluted_Psi]
display(plot(x_vec, real.(module_Psi), label=labels, xlabel=L"x", ylabel=L"ψ(x)", title="First 3 eigenfunctions "*L"ψ_1 \, ψ_2 \, ψ_3"))


# Time evolution from δ(x' - 0.5)
psi_0 = initialize_from_delta(N, dx)
println("Norm of Ψ_0 as delta: ", L2_integral(psi_0, dx))
alpha_n = [integrate_product(psi_exact[:,i], psi_0, N, dx) for i in 2:N]
println("α_n: ",  alpha_n)
@show evoluted_Psi = time_evolution(alpha_n, psi_exact, lambda_exact, 1, N)
println("Norm of Ψ(t)", L2_integral(evoluted_Psi, dx))
# I have no idea what happens for t>0
if false
    anim = @animate for k in 1:200
        evoluted_Psi = time_evolution(alpha_n, psi_exact, lambda_exact, k/100000, N) #t_bar*(1-0.01*(k-50))
        module_Psi = [conj(psi_i)*psi_i for psi_i in evoluted_Psi]
        plot(x_vec, real.(module_Psi), label=L"|Ψ|^2", ylim=(-0.05, 10.1), xlabel=L"x", ylabel=L"|Ψ|^2(x)", title="Evolution of "*L"δ(x-0.5)"*@sprintf " at time %.00005f" k/100000)
    end
    display(gif(anim, "/home/frossi/ComputationalPhysics/Assignment_2/Time_evolution_delta.gif", fps=60))
end





##############################
# BOX WITH POTENTIAL BARRIER #
##############################

function solution_FDM_box(V_n, N, dx)
    H = SymTridiagonal(2*ones(N), -ones(N-1)) / (dx^2)
    H = H + diagm(V_n)
    # display(H)
    lambda_vec = eigvals(H)
    psi_vec = eigvecs(H)
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

function initialize_from_psi1(psi_barrier)
    # This is eigenfunction ψ_1 
    psi_0 = psi_barrier[:,2]
    return psi_0
end

function t_evolution_euler(psi_0, N, dx, V_n, dt)
    H = SymTridiagonal(2*ones(N), -ones(N-1)) / (dx^2)
    H = H + diagm(V_n)
    H = diagm(ones(N)) - 1im * dt * H
    psi_0_evol = zeros(N+2)im
    psi_0_evol[2:end-1] = H*psi_0[2:end-1]
    return psi_0_evol
end

function t_evolution_CN(psi_0, N, dx, V_n, dt)
    H = SymTridiagonal(2*ones(N), -ones(N-1)) / (dx^2)
    H = H + diagm(V_n)
    H = inv(diagm(ones(N)) + 0.5im * dt * H) * (diagm(ones(N)) - 0.5im * dt * H)
    psi_0_evol = zeros(N+2)im
    psi_0_evol[2:end-1] = H*psi_0[2:end-1]
    return psi_0_evol
end
    

# Checking the results
V_n = zeros(N-1)
lambda_barrier, psi_barrier = solution_FDM_box(V_n, N-1, dx)

v0 = 1e3
V_n = potential_barrier(v0, N, dx)
lambda_barrier, psi_barrier = solution_FDM_box(V_n, N-1, dx)
println("First 10 λ_n from box with barrier of v_0 = $v0 : ", lambda_barrier[1:10])
println("Last 10 λ_n from box with barrier of v_0 = $v0 : ", lambda_barrier[90:99])
plot(1:20, lambda_exact[1:20], xlabel=L"n", ylabel=L"λ_n", label=L"λ_n, v_0 = 0"*", exact", marker = :circle, title="First 20 eigenvalues depending on "*L"n")
plot!(1:20, lambda_vec[1:20], label=L"λ_n, v_0 = 0", marker = :circle)
display(plot!(1:20, lambda_barrier[1:20], label=L"λ_n, v_0 = 10^3", marker = :circle))
display(plot(2:N-1, (lambda_barrier[2:end].-lambda_barrier[1:end-1])./lambda_barrier[2:end], marker = :circle,))

# Plotting first 3 eigenfunctions
plotting=[2,3,8]
labels=["Numerical "*L"ψ_1" "Numerical "*L"ψ_2" "Numerical "*L"ψ_7"]
display(plot(x_vec, psi_barrier[:,plotting], xlabel=L"x", ylabel=L"ψ(x)", label=labels, title="Selection of eigenfunctions with "*L"v_0 = 10^3"))
# println(lambda_barrier)
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

#### Root finding ###
root_f = [root_function(λ_i, v0) for λ_i in 1:v0-1]
display(plot( 1:v0-1, root_f, xlabel="λ", label="f(λ)", title="Root function"))

# Here the root_finding algorithm, if I had one
#2/3 x cos(x/3) (6 sin(x/3) + x cos(x/3)) = 0
 
#Task 3.6 if we lower v0, we get less states with λ < v0. When cos = 0 we get a sinh which is always positive



### Step by step evolution ###
println("Euler method")
for i in 1:100
    psi_0 = initialize_from_psi1(psi_barrier)
    global psi_0 = t_evolution_euler(psi_0, N-1, dx, V_n, 0.0001*i)
    println("At time : ",0.0001*i, " squared norm is : ", L2_integral(psi_0, dx))
end
#NOTE L2 norm is increasing

println("Crank nicolson method")
for i in 1:10
    psi_0 = initialize_from_psi1(psi_barrier)
    global psi_0 = t_evolution_CN(psi_0, N-1, dx, V_n, 0.0001*i)
    println("At time : ",0.0001*i, " squared norm is : ", L2_integral(psi_0, dx))
end
if false
    anim = @animate for i in 1:10
        psi_0 = initialize_from_psi1(psi_barrier)
        global psi_0 = t_evolution_CN(psi_0, N-1, dx, V_n, 100*i)
        module_Psi = [conj(psi_i)*psi_i for psi_i in psi_0]
        plot(x_vec, real.(module_Psi), label="Squared module", xlabel="Index", ylabel="|Ψ|^2", title="Squared module at time *"*string(0.01*i))
    end
    display(gif(anim, "/home/frossi/ComputationalPhysics/Assignment_2/Time_evolution_barrier_CN.gif", fps=60))
end
# If I push it to large times it won't change.
# Not even when changing psi_0 to superposition.
# Not even when making Psi_0 Im
# Not even with the delta


#####################
# PERIODIC DETUNING #
#####################

function double_potential_barrier(V0, V1, N, dx)
    V_n = zeros(N-1)
    for i in 1:N-1
        if 1/3 <= i*dx < 2/3
            V_n[i] = V0
        elseif 2/3 <= i*dx < 1
            V_n[i] = V1
        end
    end
    return V_n
end

function evaluate_τ(V_n, ψ, N, dx)
    H = SymTridiagonal(2*ones(N), -ones(N-1)) / (dx^2)
    H = H + diagm(V_n)
    ψ_test = zeros(N+2)
    ψ_test[2:end-1] = H*ψ[2:end-1,3]

    return integrate_product(ψ[:,2], ψ_test, N, dx)
end


λ1_list = []
λ2_list = []
v_list = []

plotting = [2,3]
labels = ["ψ1" "ψ2"]

for v1 in -500:500
    local V_n = double_potential_barrier(v0, v1, N, dx)
    lambda_double, psi_double = solution_FDM_box(V_n, N-1, dx)
    push!(v_list,  v1)
    push!(λ1_list, lambda_double[1])
    push!(λ2_list, lambda_double[2])
    if v1 == -500
        display(plot(x_vec, psi_double[:,plotting], label=labels, title="First eigenfunctions at V1 = -500"))
    elseif v1 == 500
        display(plot(x_vec, psi_double[:,plotting], label=labels, title="First eigenfunctions at V1 = 500"))
    elseif v1 == 0
        println("ΔE at V1=0 : ", lambda_double[2]-lambda_double[1])
    end
end

plot(v_list, λ1_list, xlabel="V1", label=" λ0", title="Eigenvalues depending on V1")
display(plot!(v_list, λ2_list, label="λ1"))

display(plot(v_list, λ2_list .- λ1_list, xlabel="V1", label="λ1-λ0", title="Eigenvalues difference depending on V1"))

# Finding τ
V_n = double_potential_barrier(v0, 0.0, N, dx)
lambda_ref, psi_ref = solution_FDM_box(V_n, N-1, dx)
println("τ with v1 = 0 : ", evaluate_τ(V_n, psi_ref, N-1, dx))

τ_list = []
v_list = []
for v1 in -500:500
    local V_n = double_potential_barrier(v0, v1, N, dx)
    push!(τ_list, evaluate_τ(V_n, psi_ref, N-1, dx))
    push!(v_list,  v1)
end
display(plot(v_list, τ_list, xlabel="V1", label="τ(V1)", title="τ depending on V1"))
