using LinearAlgebra
using Plots
using SpecialFunctions
using LaTeXStrings
using Printf
using Arpack
using Random
using LinearSolve
using IterativeSolvers

function build_T(N)
    # Builds the transformation matrix T
    T_mat = zeros(N,N)
    for i in 1:N-2
        T_mat[i,i+1] = 1.0
        T_mat[i,i+2] = 1.0
    end
    T_mat[N-1,N] = 1.0
    T_mat[1:2,N] .= 1.0
    T_mat[1,N-1] = 1.0

    return 0.25*(T_mat+ T_mat')
end

function build_disjoint_T(N1, N2)
    T_mat = zeros(N1+N2, N1+N2)
    T_mat[1:N1, 1:N1] = build_T(N1)
    T_mat[N1+1:end, N1+1:end] = build_T(N2)
    return T_mat
end

function plot_bars(values, title)
    N = length(values)
    x = 1:N
    bar(x, values, ylim=(0,1), xlabel=L"n", ylabel=L"V(n)", title=title, legend=false)
end

function time_evolution_n(V, T_mat, n_steps, name)
    anim = @animate for k in 0:n_steps
        if k>0
            V = T_mat * V
        end
        plot_bars(V,L"t="*string(k)*L"\, dt")    
    end
    display(gif(anim, "/home/frossi/ComputationalPhysics/Exam/Prob2/$name.gif", fps=4))
    return V
end

function iterative_extreme_eigvals(A, n_extremes)
    # Arpack uses Lanczos algorithm for symmetric matrices
    eigvals, _ = Arpack.eigs(A, nev=n_extremes, which=:BE)
    return eigvals
end

function solve_linear_problem(A, b, method)
    if method=="invert"
        return inv(A)*b
    elseif method=="LU"
        prob = LinearProblem(A, b)
        sol = solve(prob, LUFactorization()) 
        return sol.u
    elseif method=="CG"
        return cg(A, b)
    elseif method=="Krylov_GMRES"
        prob = LinearProblem(A, b)
        sol = solve(prob, KrylovJL_GMRES()) 
        return sol.u
    end
end

N=21
x = 1:N
### Building matrix T [Task 2a]
T_mat = build_T(N)
display(T_mat)

### Running time evolution [Task 2b]
V = zeros(N)
V[1] = 1.0
bar(x, V, ylim=(0,1), label=L"t=0", xlabel=L"n", ylabel=L"V(n)", legend=true)
V = T_mat * V
bar!(x, V, label=L"t = dt")
V = T_mat * V
bar!(x, V, label=L"t=2\,dt")
V = T_mat^98 * V
display(bar!(x, V, label=L"t=100\,dt"))

V = zeros(N)
V[1] = 1.0
V = time_evolution_n(V, T_mat, 30, "time_evol_100000")
display(plot_bars(V, "Time evolution, V(t=0)=1,0,0,0..."))

### Smaller and bigger eigvals [Task 2c]
println("Iterative extreme eigvals: ", iterative_extreme_eigvals(T_mat, 2))
println("Full list of eigvals:")
@show evals = eigvals(T_mat)
evecs = eigvecs(T_mat)
# ANSWER: The eigvals are in pairs and then a last one isolated at 1, only one


### Comparing at t=100 dt [Task 2d]
# V = 1, 0, 0, 0, ...
if true     # PUT TRUE TO RUN
    V = zeros(N)
    V[1] = 1.0
    V = time_evolution_n(V, T_mat, 100, "time_evol_100000")
    display(plot_bars(V, ""))

    # V = some points...
    V = zeros(N)
    for i in [2,7,8,15,18]
        V[i] = 1.0/5
    end
    V = time_evolution_n(V, T_mat, 100, "time_evol_somepoints")
    display(plot_bars(V, ""))

    # V random
    V = rand(N)
    V = V / sum(V)
    V = time_evolution_n(V, T_mat, 100, "time_evol_random")
    display(plot_bars(V, ""))
end  
# ANSWER: They are related to the last eigenvector.
#         This is because V(100dt) = T^100 V(0) ≈ |v_N><v_N|V(0)>


### Disjoint network [Task 2e]
disj_T = build_disjoint_T(11,10)
println("Iterative extreme eigvals: ", iterative_extreme_eigvals(disj_T, 4))
evals = eigvals(disj_T)
display(evals)
evecs = eigvecs(disj_T)
# ANSWER: Two sets of eigvals and eigvecs, each set in blocks
#         Eigvals in pairs. 11 in 5 pairs and one "1".
#                           10 in 4 pairs and one "1" and one "0"

if true     # PUT TRUE TO RUN
    V = zeros(N)
    V[1] = 1.0
    V = time_evolution_n(V, disj_T, 30, "time_evol_disj_100000")
    display(plot_bars(V, ""))
    # This goes to 1/11 in the first block

    # V = 0,..0, 1_12,0,..0
    V = zeros(N)
    V[12] = 1.0
    V = time_evolution_n(V, disj_T, 30, "time_evol_disj_000100")
    display(plot_bars(V, ""))
    # This goes to 1/10 in the second block

    # V 1, 0, ..0, 1_12, 0, .., 0
    V = zeros(N)
    V[1] = 0.5
    V[12] = 0.5
    V = time_evolution_n(V, disj_T, 30, "time_evol_disj_100100")
    display(plot_bars(V, ""))

    # V 0.75, 0, ..0, 0.25_12, 0, .., 0
    V = zeros(N)
    V[1] = 0.75
    V[12] = 0.25
    V = time_evolution_n(V, disj_T, 30, "time_evol_disj_300100")
    display(plot_bars(V, ""))
end

### Trace back V[1/i] [Task 2f]
V = [1.0/i for i in 1:N]
methods = ["invert", "LU", "CG", "Krylov_GMRES"]
for met in methods
    if met == "invert"
        global V_exact = solve_linear_problem(T_mat, V, met)
        println(V_exact)
    end
    println(met)
    new_V = solve_linear_problem(T_mat, V, met) # precompile
    @time solve_linear_problem(T_mat, V, met)
    println("Error is : ", norm(new_V - V_exact))
end

# Now t-5dt
for met in methods
    V = [1.0/i for i in 1:N]
    println(met)
    solve_linear_problem(T_mat, V, met) # to precompile
    @time begin
        for _ in 1:5
            V = solve_linear_problem(T_mat, V, met)
        end
    end
    if met == "invert"
        global V_exact = V
        println(V)
    end
    println("Error is : ", norm(V - V_exact))
end