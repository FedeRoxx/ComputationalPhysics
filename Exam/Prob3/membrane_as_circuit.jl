using LinearAlgebra
using Plots
using SpecialFunctions
using LaTeXStrings
using Printf

function initialize_gaussian(x, x0, σ)
    gaussian = exp.(-0.5 * ((x .- x0) / σ).^2) / (σ * sqrt(2 * π))
    return gaussian / sum(gaussian)
end

function trapezoidal_integral(f, dx)
    return (sum(f)-0.5*(f[1]-f[end])) * dx
end

function euler_explicit_step(v, α, β)
    N = length(v)
    new_v = (1-β)*v
    for i in 2:N-1
        new_v[i] += α*(v[i+1]+v[i-1]-2*v[i])
    end
    # Enforce Neumann boundary conditions
    new_v[1] += 2*α*(v[2]-v[1])
    new_v[N] += 2*α*(v[N-1]-v[N])
    return new_v
end

function run_euler_explicit(v, α, β, n_steps)
    v_list = []
    anim = @animate for k in 1:n_steps
        v = euler_explicit_step(v, α, β)
        plot_potential(x, v, "Euler explicit, step "*string(k))
        push!(v_list, v)
    end
    display(gif(anim, "/home/frossi/ComputationalPhysics/Exam/Prob3/Euler_explicit.gif", fps=10))
    display(plot_potential(x, v, "Euler explicit"))
    return v_list
end

function plot_potential(x, v, title)
    plot(x, v, ylim=(-0.005, 5), xlabel=L"x", ylabel=L"V(x)", title=title, legend=false)
end

function analytical_unbound(x, x0, λ, τ, V0, t)
    C = 4*λ^2/τ * t
    v_exact = exp.(-(x .- x0).^2 / C .- t/τ)    
    return v_exact * V0 / sqrt(π*C)
end

function run_exact_unbound(x, x0, λ, τ, V0, dt, n_steps)
    v_list = []
    anim = @animate for k in 1:n_steps
        v = analytical_unbound(x, x0, λ, τ, V0, k*dt)
        plot_potential(x, v, "Analytical unbound, step "*string(k))
        push!(v_list, v)
    end
    display(gif(anim, "/home/frossi/ComputationalPhysics/Exam/Prob3/Exact.gif", fps=10))
    display(plot_potential(x, v_list[end], "Analytical unbound"))
    return v_list
end

function euler_implicit_step(A, v)
    return A \ v
end

function run_euler_implicit(v, α, β, n_steps)
    N = length(v)

    diag = ones(N) * (1 + 2α + β)
    upp_diag = -α* ones(N-1)
    upp_diag[1] = -2α
    low_diag = -α* ones(N-1)
    low_diag[end] = -2α

    A = Tridiagonal(low_diag, diag, upp_diag)
    v_list = []
    anim = @animate for k in 1:n_steps
        v = euler_implicit_step(A, v)
        plot_potential(x, v, "Euler implicit, step "*string(k))
        push!(v_list, v)
    end
    display(gif(anim, "/home/frossi/ComputationalPhysics/Exam/Prob3/Euler_implicit.gif", fps=10))
    display(plot_potential(x, v, "Euler implicit"))
    return v_list
end

function run_CN(v, α, β, n_steps)
    N = length(v)

    diag_plus = ones(N) * (1 + α + β/2)
    diag_minus = ones(N) * (1 - α - β/2)
    upp_diag = -α/2* ones(N-1)
    upp_diag[1] = -α
    low_diag = -α/2* ones(N-1)
    low_diag[end] = -α
    A = Tridiagonal(low_diag, diag_plus, upp_diag)
    B = Tridiagonal(-low_diag, diag_minus, -upp_diag)

    C = inv(A)*B  # This is only evaluated once

    v_list = []
    anim = @animate for k in 1:n_steps
        v = C*v
        plot_potential(x, v, "Crank Nicolson, step "*string(k))
        push!(v_list, v)
    end
    display(gif(anim, "/home/frossi/ComputationalPhysics/Exam/Prob3/Crank_Nicolson.gif", fps=10))
    display(plot_potential(x, v, "Crank Nicolson"))
    return v_list
end

function plot_list(v_lists, labels, k)
    p = plot(ylim=(-0.005, 5), label=labels[1], xlabel=L"x", ylabel=L"V(x)", title="All methods compared, step "*string(k), legend=true)
    for (i, v_list) in enumerate(v_lists)
        if labels[i] == "Exact unbound"
            plot!(p, x, v_list[k], label=labels[i])
        else
            scatter!(p, x, v_list[k], label=labels[i], marker=:x)
        end
    end
    return plot!(p)
end

function compare_with_gif(v_lists, labels)
    n_steps = length(v_lists[1])
    anim = @animate for k in 1:n_steps # 
        plot_list(v_lists, labels, k)
    end
    display(gif(anim, "/home/frossi/ComputationalPhysics/Exam/Prob3/Comparison.gif", fps=10))
end




dt = 0.00015

L = 1.0  # Interval is [0,L]
N = 51 # n of points
@show dx = L/(N-1)
x = [i*dx for i in 0:N-1]
x0 = L/2
λ = 1.0
τ = 1.0
β = dt / τ
@show α = λ^2*β/dx^2

n_steps = 500

v = initialize_gaussian(x, x0, 0.01)
@show V0 = trapezoidal_integral(v, dx)
v = v / V0
display(plot_potential(x, v, "Before"))

v_list_EE = run_euler_explicit(v, α, β, n_steps)

v_list_EI = run_euler_implicit(v, α, β, n_steps)

v_list_CN = run_CN(v, α, β, n_steps)

v_list_exact = run_exact_unbound(x, x0, λ, τ, 1.0, dt, n_steps)
# display(plot_potential(x, v_list_exact[end], "After EXACT"))

compare_with_gif([v_list_EE, v_list_EI, v_list_CN, v_list_exact], ["Explicit euler", "Implicit euler", "Crank Nicolson", "Exact unbound"])