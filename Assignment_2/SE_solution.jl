using LinearAlgebra
using Plots
using SpecialFunctions
using LaTeXStrings
using Printf
import Pkg; Pkg.add("Optim")
using Optim

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
println("L2 norm of psi_exact_1 - psi_numerical_1 : ", min(L2_integral(psi_exact[:,2]-psi_vec[:,2], dx), L2_integral(psi_exact[:,2]+psi_vec[:,2], dx)))
println("L2 norm of psi_exact_2 - psi_numerical_2 : ", min(L2_integral(psi_exact[:,3]-psi_vec[:,3], dx), L2_integral(psi_exact[:,3]+psi_vec[:,3], dx)))

# Plotting eigenval in function of n
plot(n_vec, lambda_exact, xlabel=L"n", ylabel=L"λ_n", label=L"λ_n"*" exact", title="Comparison of eigenvalues as function of "*L"n")
display(plot!(n_vec, lambda_vec, label=L"λ_n"*" numerical"))

# Relative difference
display(plot(n_vec, (lambda_exact.-lambda_vec)./lambda_exact, label=L"λ_n"*" relative difference"))

# Plotting numerical vs exact_solution eigenvalues
display(plot(lambda_exact, lambda_vec, label="Numerical"))





################################
# EXPANSIONS IN EIGENFUNCTIONS #
################################

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

# Animation
if false
    anim = @animate for k in 0:200
        evoluted_Psi = time_evolution(alpha_n, psi_exact, lambda_exact, 1/(pi)+0.00001*k, N) #t_bar*(1-0.01*(k-50))
        module_Psi = [conj(psi_i)*psi_i for psi_i in evoluted_Psi]
        plot(x_vec, real.(module_Psi), label=L"|Ψ|^2", ylim=(-0.05, 10.1), xlabel=L"x", ylabel=L"|Ψ|^2(x)", title="Evolution of "*L"δ(x-0.5)"*@sprintf " at time 1/π + %.5f" 0.00001*k)
    end
    display(gif(anim, "/home/frossi/ComputationalPhysics/Assignment_2/Time_evolution_delta.gif", fps=30))
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
    # This is eigenfunction ψ_1+ψ_2 / sqrt(2)
    psi_0 = (psi_barrier[:,2]+psi_barrier[:,3])/sqrt(2)
    return psi_0
end

function root_function(λ)
    # if λ > v0
    #     throw("Root function is defined for 0 < λ , v0")
    # end
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

function find_nearest_minimum(f, λ_i)
    result = optimize(x -> f(x[1]), [λ_i], NelderMead(; parameters = Optim.FixedParameters(),
    initial_simplex = Optim.AffineSimplexer()), Optim.Options(f_tol = 1e-6, g_tol = 1e-10))
    min_point = Optim.minimizer(result)
    # println(Optim.converged(result), Optim.iterations(result))
    return min_point[1]
end

function find_root_from_min(f, min_λ, α)
    x_minus = min_λ
    x_plus = α* min_λ
    if x_plus < 0
        return find_root_from_min(f, min_λ, α)
    end
    tol = 1e-6
    while abs(x_plus - x_minus) / 2 > tol
        x_mean = (x_plus + x_minus) / 2
        if abs(f(x_mean)) < tol
            return x_mean
        elseif f(x_mean) > 0
            x_plus = x_mean
        else
            x_minus = x_mean
        end
    end
    return (x_plus + x_minus) / 2
end

function find_roots(lambda_barrier, v0)
    λ_list = [lambda_barrier[i] for i in 1:2:length(lambda_barrier)-1 if lambda_barrier[i] < v0]
    min_list = [find_nearest_minimum(root_function, λ_i) for λ_i in λ_list]
    roots = []
    for i in eachindex(λ_list)
        push!(roots, find_root_from_min(root_function, min_list[i], 0.99))
        push!(roots, find_root_from_min(root_function, min_list[i], 1.01))
    end
    return roots
end


# Checking with v_0 = 0
V_n = zeros(N-1)
lambda_barrier, psi_barrier = solution_FDM_box(V_n, N-1, dx)

# Now v_0 = 1000
N = 100
dx = 1/N
v0 = 1e3
V_n = potential_barrier(v0, N, dx)
lambda_barrier, psi_barrier = solution_FDM_box(V_n, N-1, dx)
println("First 10 λ_n from box with barrier of v_0 = $v0 : ", lambda_barrier[1:10])
plot(1:20, lambda_exact[1:20], xlabel=L"n", ylabel=L"λ_n", label=L"λ_n, v_0 = 0"*", exact", marker = :circle, title="First 20 eigenvalues depending on "*L"n")
plot!(1:20, lambda_vec[1:20], label=L"λ_n, v_0 = 0", marker = :circle)
display(plot!(1:20, lambda_barrier[1:20], label=L"λ_n, v_0 = 10^3", marker = :circle))
display(plot(2:N-1, (lambda_barrier[2:end].-lambda_barrier[1:end-1])./lambda_barrier[2:end], marker = :circle,))

# Plotting first 3 eigenfunctions
plotting=[2,3,8]
labels=["Numerical "*L"ψ_1" "Numerical "*L"ψ_2" "Numerical "*L"ψ_7"]
display(plot(x_vec, psi_barrier[:,plotting], xlabel=L"x", ylabel=L"ψ(x)", label=labels, title="Selection of eigenfunctions with "*L"v_0 = 10^3"))

# Time evolution from superposition
psi_0 = initialize_from_superposition(psi_barrier)
display(plot(x_vec, psi_0, title="Initial condition"))
println("Norm of Ψ_0 as (ψ_1 + ψ_2) / √2: ", L2_integral(psi_0, dx))
alpha_n = [integrate_product(psi_barrier[:,i], psi_0, N, dx) for i in 2:N]
println("α_n: ",  alpha_n)
evoluted_Psi = time_evolution(alpha_n, psi_barrier, lambda_barrier, 1, N)
println("Norm of Ψ(t)", L2_integral(evoluted_Psi, dx))

t_bar=π/(lambda_barrier[2]-lambda_barrier[1])
println("τ for evolution of superposition initial condition: ", t_bar)
if false
anim = @animate for k in 1:200
    evoluted_Psi = time_evolution(alpha_n, psi_barrier, lambda_barrier, t_bar*k/100, N) #t_bar*(1-0.01*(k-50))
    module_Psi = [conj(psi_i)*psi_i for psi_i in evoluted_Psi]
    plot(x_vec, real.(module_Psi), ylim=(-0.05, 6), label=L"|Ψ|^2", xlabel=L"x", ylabel=L"|Ψ|^2(x)", title=@sprintf "Squared module at time %.2f τ" k/100 )
end
display(gif(anim, "/home/frossi/ComputationalPhysics/Assignment_2/Time_evolution_barrier.gif", fps=60))
end

#### Root finding ###
root_f = [root_function(λ_i) for λ_i in 1:v0]
display(plot( 1:v0, root_f, xlabel=L"λ", label=L"f(λ)", title="Root function at "*L"v_0 = 10^3"))

v0=1e5
root_f = [root_function(λ_i) for λ_i in 1:v0]
display(plot( 1:v0, root_f, xlim=(0,2.5e4), xlabel=L"λ", label=L"f(λ), v_0 = 10^5", title="Root function at "*L"v_0 = 10^5"))

v0 = 24
root_f = [root_function(λ_i/10) for λ_i in 1:10*v0]
x = [λ_i/10 for λ_i in 1:10*v0]
plot(x, root_f, xlim=(0,25), xticks = 0:1:25, xlabel=L"λ", label=L"f(λ), v_0 = 24", title="Root function around "*L"v_0 = 22")
v0 = 23
root_f = [root_function(λ_i/10) for λ_i in 1:10*v0]
x = [λ_i/10 for λ_i in 1:10*v0]
plot!(x, root_f, xlim=(0,25), label=L"f(λ), v_0 = 23")
v0 = 22
root_f = [root_function(λ_i/10) for λ_i in 1:10*v0]
x = [λ_i/10 for λ_i in 1:10*v0]
display(plot!( x, root_f, xlim=(0,25), label=L"f(λ), v_0 = 22"))

v0 = 1e3

# Here the root_finding algorithm
println("Numerical roots: ", find_roots(lambda_barrier, v0))


### Step by step evolution ###
println("Euler method")
#Without animation
psi_0 = initialize_from_psi1(psi_barrier)
for i in 1:10
    global psi_0 = t_evolution_euler(psi_0, N-1, dx, V_n, 0.00003)
    println("At time : ",0.00003*i, " squared norm is : ", L2_integral(psi_0, dx))
end
#Animation
if false
    psi_0 = initialize_from_psi1(psi_barrier)
    anim = @animate for i in 1:100
        # psi_0 = initialize_from_psi1(psi_barrier)
        global psi_0 = t_evolution_euler(psi_0, N-1, dx, V_n, 0.00003)
        module_Psi = [conj(psi_i)*psi_i for psi_i in psi_0]
        plot(x_vec, real.(module_Psi), ylim=(-0.05,5), label=L"|Ψ|^2", xlabel=L"x", ylabel=L"|Ψ|^2(x)", title=@sprintf "Squared module, t= %.5f L2 norm: %.2f " 0.00003*i L2_integral(psi_0, dx))
    end
    display(gif(anim, "/home/frossi/ComputationalPhysics/Assignment_2/Time_evolution_barrier_FE.gif", fps=60))
end


N=100
dx = 1/N
x_vec = dx*[i for i in 0:N]
V_n = potential_barrier(v0, N, dx)


println("Crank nicolson method")
psi_0 = initialize_from_delta(N, dx)
# Wihtout animation
for i in 1:10
    global psi_0 = t_evolution_CN(psi_0, N-1, dx, V_n, 0.0001)
    println("At time : ",0.0001*i, " squared norm is : ", L2_integral(psi_0, dx))
end
#Animation from delta
if false
    psi_0 = initialize_from_delta(N, dx)
    anim = @animate for i in 1:500
        # psi_0 = initialize_from_delta(N, dx)
        global psi_0 = t_evolution_CN(psi_0, N-1, dx, V_n, 0.00001)
        module_Psi = [conj(psi_i)*psi_i for psi_i in psi_0]
        plot(x_vec, real.(module_Psi), ylim=(-0.05,15),label=L"|Ψ|^2", xlabel=L"x", ylabel=L"|Ψ|^2(x)", title=@sprintf "Squared module, t= %.4f L2 norm: %.2f " 0.00001*i L2_integral(psi_0, dx))
    end
    display(gif(anim, "/home/frossi/ComputationalPhysics/Assignment_2/Time_evolution_barrier_CN_delta.gif", fps=60))
end
# Animation from ψ_1
if false
    psi_0 = initialize_from_psi1(psi_barrier)
    anim = @animate for i in 1:1000
        # psi_0 = initialize_from_psi1(psi_barrier)
        global psi_0 = t_evolution_CN(psi_0, N-1, dx, V_n, 0.001)
        module_Psi = [conj(psi_i)*psi_i for psi_i in psi_0]
        plot(x_vec, real.(module_Psi), ylim=(-0.05,5),label=L"|Ψ|^2", xlabel=L"x", ylabel=L"|Ψ|^2(x)", title=@sprintf "Squared module, t= %.3f L2 norm: %.2f " 0.001*i L2_integral(psi_0, dx))
    end
    display(gif(anim, "/home/frossi/ComputationalPhysics/Assignment_2/Time_evolution_barrier_CN.gif", fps=60))
end

N=100
dx = 1/N
x_vec = dx*[i for i in 0:N]
V_n = potential_barrier(v0, N, dx)

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

function build_N(ε0, t, ω, α_τ)
    Nk = zeros(2,2)im
    Nk[1,2] = exp(-1im*t*ε0)
    Nk[2,1] = exp(1im*t*ε0)
    return Nk * α_τ *ε0*sin(ω*t)
end

function solve_volterra(ψ_0, max_k, dt, ε0, ω, α_τ)
    f_list = [ψ_0]
    N_list = []
    for k in 1:max_k
        Nk = -1im*dt*build_N(ε0, k*dt, ω, α_τ)
        Mk = inv(diagm(ones(2))-0.5*Nk)
        f0 = copy(ψ_0)
        for l in eachindex(N_list)
            f0 += N_list[l]*f_list[l+1]
        end
        push!(f_list, Mk*f0)
        push!(N_list, Nk)
    end
    return f_list
end

# Double barrier at different v_1
λ1_list = []
λ2_list = []
v_list = []

plotting = [2,3]
labels = [L"ψ_1" L"ψ_2"]

for v1 in -500:500
    local V_n = double_potential_barrier(v0, v1, N, dx)
    lambda_double, psi_double = solution_FDM_box(V_n, N-1, dx)
    push!(v_list,  v1)
    push!(λ1_list, lambda_double[1])
    push!(λ2_list, lambda_double[2])
    if v1 == -500
        display(plot(x_vec, psi_double[:,plotting], xlabel=L"x", ylabel=L"ψ(x)",label=labels, title="First eigenfunctions at "*L"v_1 = -500"))
    elseif v1 == 500
        display(plot(x_vec, psi_double[:,plotting], xlabel=L"x", ylabel=L"ψ(x)", label=labels, title="First eigenfunctions at "*L"v_1 = 500"))
    elseif v1 == 0
        global ε0 = lambda_double[2]-lambda_double[1]
        println("ΔE at V_1 = 0 : ", ε0)
    end
end

plot(v_list, λ1_list, xlabel=L"v_1", label=L" λ_0", title="Eigenvalues depending on "*L"v_1")
display(plot!(v_list, λ2_list, label=L"λ_1"))

display(plot(v_list, λ2_list .- λ1_list, xlabel=L"v_1", label=L"λ_1-λ_0", title="Eigenvalues difference depending on "*L"v_1"))

# Checks for tau when v1=0 
V_n = double_potential_barrier(v0, 0.0, N, dx)
lambda_ref, psi_ref = solution_FDM_box(V_n, N-1, dx)
display(plot(x_vec, psi_ref[:,2:3]))
println("τ with v1 = 0 : ", evaluate_τ(V_n, psi_ref, N-1, dx))

# Tau varying v1
τ_list = []
v_list = []
for v1 in -500:500
    local V_n = double_potential_barrier(v0, v1, N, dx)
    push!(τ_list, evaluate_τ(V_n, psi_ref, N-1, dx))
    push!(v_list,  v1)
end
# Finding m in τ = m*v1
# println(τ_list ./ v_list)

display(plot(v_list, τ_list, xlabel=L"v_1", label=L"τ(v_1)", title=L"τ "*" as a function of "*L" v_1"))

# Probability around T1
dt = 5
max_k = 1200
ψ_n = solve_volterra([1.0+0.0im,0.0im], max_k, dt, ε0, ε0, 0.02)
t_vec = [dt*i for i in 0:max_k]
# display(ψ_n)
prob = [abs(ψ_n[k][2])^2 for k in 1:max_k+1]
prob_exact = [sin(0.02*ε0*i*dt/2)^2 for i in 0:max_k]
plot(t_vec, prob_exact , label=L"p(t) = "*"sin"*L"^2  \left( \frac{t\tau}{2}\right)", xlabel=L"t", ylabel=L"p(t)", title=L"p(t) "*" for "*L"t "*" ~ "*L" T_1" )
display(plot!(t_vec, prob , label="Numerical "*L"p(t)"))

# Probability around T2
dt = 500
max_k = 500
ψ_n = solve_volterra([1.0+0.0im,0.0im], max_k, dt, ε0, ε0, 0.02)
t_vec = [dt*i for i in 0:max_k]
# display(ψ_n)
prob = [abs(ψ_n[k][2])^2 for k in 1:max_k+1]
prob_exact = [sin(0.02*ε0*i*dt/2)^2 for i in 0:max_k]
plot(t_vec, prob_exact, label=L"p(t) = "*"sin"*L"^2  \left( \frac{t\tau}{2}\right)", xlabel=L"t", ylabel=L"p(t)", title=L"p(t) "*" for "*L"t "*" ~ "*L" T_2")
# plot!([14300],[0.0], marker = :circle)
display(plot!(t_vec, prob ,label="Numerical "*L"p(t)"))

# Frequency out of resonance
dt = 5
max_k = 1200
ω = ε0 *0.95
ψ_n = solve_volterra([1.0+0.0im,0.0im], max_k, dt, ε0, ω, 0.02)
t_vec = [dt*i for i in 0:max_k]
# display(ψ_n)
prob = [abs(ψ_n[k][2])^2 for k in 1:max_k+1]
prob_exact = [sin(0.02*ε0*i*dt/2)^2 for i in 0:max_k]
plot(t_vec, prob_exact , label=L"p(t) = "*"sin"*L"^2  \left( \frac{t\tau}{2}\right)", xlabel=L"t", ylabel=L"p(t)", title=L"p(t) "*" for "*L"t "*" ~ "*L" T_1 \quad (ω=0.95 λ_0)" )
display(plot!(t_vec, prob , label="Numerical "*L"p(t)"))

dt = 500
max_k = 500
ψ_n = solve_volterra([1.0+0.0im,0.0im], max_k, dt, ε0, ω, 0.02)
t_vec = [dt*i for i in 0:max_k]
# display(ψ_n)
prob = [abs(ψ_n[k][2])^2 for k in 1:max_k+1]
prob_exact = [sin(0.02*ε0*i*dt/2)^2 for i in 0:max_k]
plot(t_vec, prob_exact, label=L"p(t) = "*"sin"*L"^2  \left( \frac{t\tau}{2}\right)", xlabel=L"t", ylabel=L"p(t)", title=L"p(t) "*" for "*L"t "*" ~ "*L" T_2 \quad (ω=0.95 λ_0)")
# plot!([14300],[0.0], marker = :circle)
display(plot!(t_vec, prob ,label="Numerical "*L"p(t)"))




# Resonance but large τ
α_τ = 0.75
dt = 5
max_k = 2500
ψ_n = solve_volterra([1.0+0.0im,0.0im], max_k, dt, ε0, ε0, α_τ)
t_vec = [dt*i for i in 0:max_k]
# display(ψ_n)
prob = [abs(ψ_n[k][2])^2 for k in 1:max_k+1]
prob_exact = [sin(α_τ*ε0*i*dt/2)^2 for i in 0:max_k]
plot(t_vec, prob_exact , label=L"p(t) = "*"sin"*L"^2  \left( \frac{t\tau}{2}\right)", xlabel=L"t", ylabel=L"p(t)", title=L"p(t)   "*L"\quad(τ=0.8 λ_0)" )
display(plot!(t_vec, prob , label="Numerical "*L"p(t)"))

# dt = 5
# max_k = 5000
# ψ_n = solve_volterra([1.0+0.0im,0.0im], max_k, dt, ε0, ε0, α_τ)
# t_vec = [dt*i for i in 0:max_k]
# # display(ψ_n)
# prob = [abs(ψ_n[k][2])^2 for k in 1:max_k+1]
# prob_exact = [sin(α_τ*ε0*i*dt/2)^2 for i in 0:max_k]
# plot(t_vec, prob_exact, label=L"p(t) = "*"sin"*L"^2  \left( \frac{t\tau}{2}\right)", xlabel=L"t", ylabel=L"p(t)", title=L"p(t) "*" for "*L"t "*" ~ "*L" T_2 \quad (τ=0.8 λ_0)")
# # plot!([14300],[0.0], marker = :circle)
# display(plot!(t_vec, prob ,label="Numerical "*L"p(t)"))