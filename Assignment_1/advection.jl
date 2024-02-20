using LinearAlgebra
# import Pkg; Pkg.add("Plots")
using Plots
using SpecialFunctions

function next_step(C, u)
    # Looks like it should be already ok when writing the matrix C
    return C * u
end

function dirichlet(N, α, start_u, steps)
    x = [i / (N + 1) for i in 1:N]
    u_n = zeros(N)
    u_n[Int((N + 1) / 2)] = start_u
    #println(u_n)

    alpha_halved = ones(N - 1) * α / 2
    alpha_plus = ones(N) * (1 + α)
    alpha_minus = ones(N) * (1 - α)

    A = SymTridiagonal(alpha_plus, -alpha_halved)
    B = SymTridiagonal(alpha_minus, alpha_halved)
    C = inv(A) * B

    for _ in 1:steps
        u_n = C * u_n
      #  check_conservation(u_n, N)
        println(u_n)
    end
    plot(x, u_n, label="Dirichlet", show=true)
end

function neumann(N, α, start_u, steps)
    x = [i / (N + 1) for i in 1:N]
    u_n = zeros(N)
    u_n[Int((N + 1) / 2)] = start_u
    #println(u_n)

    alpha_plus = ones(N) * (1 + α)
    alpha_minus = ones(N) * (1 - α)
    alpha_halved_1 = ones(N - 1) * α / 2
    alpha_halved_1[1] = α
    alpha_halved_2 = ones(N - 1) * α / 2
    alpha_halved_2[end] = α

    A = Tridiagonal(-alpha_halved_2, alpha_plus, -alpha_halved_1)
    B = Tridiagonal(alpha_halved_2, alpha_minus, alpha_halved_1)
    C = inv(A) * B

    for _ in 1:steps
        u_n = C * u_n
        check_conservation(u_n, N)
       # println(u_n)
        # plot!(x,u_n, label="my label", show=true)
    end
    plot!(x, u_n, label="Neumann", show=true)
end

function dirichlet_space_dependent(N, α_vec, start_u, steps)
    x = [i / (N + 1) for i in 1:N]
    u_n = zeros(N)
    u_n[Int((N + 1) / 2)] = start_u
    #println(u_n)

    alpha_halved_1 = α_vec[1:end-1] / 2
    alpha_halved_2 = α_vec[2:end] / 2
    alpha_plus = ones(N) .+ α_vec
    alpha_minus = ones(N) .- α_vec

    A = Tridiagonal(-alpha_halved_2, alpha_plus, -alpha_halved_1)
    B = Tridiagonal(alpha_halved_2, alpha_minus, alpha_halved_1)
    C = inv(A) * B

    for _ in 1:steps
        u_n = C * u_n
        #check_conservation(u_n, N)
        #println(u_n)
    end
    plot(x, u_n, label="Dirichlet", show=true)
end

function neumann_space_dependent(N, α_vec, start_u, steps)
    x = [i / (N + 1) for i in 1:N]
    u_n = zeros(N)
    u_n[Int((N + 1) / 2)] = start_u
    println(u_n)

    alpha_halved_1 = α_vec[1:end-1] / 2
    alpha_halved_1[1] = α_vec[1]
    alpha_halved_2 = α_vec[2:end] / 2
    alpha_halved_2[end] = α_vec[end]
    alpha_plus = ones(N) .+ α_vec
    alpha_minus = ones(N) .- α_vec

    A = Tridiagonal(-alpha_halved_2, alpha_plus, -alpha_halved_1)
    B = Tridiagonal(alpha_halved_2, alpha_minus, alpha_halved_1)
    C = inv(A) * B

    for _ in 1:steps
        u_n = C * u_n
    #    check_conservation(u_n, N)
      #  println(u_n)
    end
    plot!(x, u_n, label="Neumann", show=true)
end

function unbounded_space_dependent(x, x0, t, D_plus, D_minus)
    #My guess is that it should be x0 = when change - when initial drop
    if x >= x0 #+ 0.2
        result = exp(-(x - x0)^2 / (4 * D_plus * t))
        return result * A_plus(t, 0.0, D_plus, D_minus) / sqrt(4 * pi * D_plus * t)
    else
        result = exp(-(x - x0)^2 / (4 * D_minus * t))
        return result * A_minus(t, 0.0, D_plus, D_minus) / sqrt(4 * pi * D_minus * t)
    end
end

function analytical_solution(x, x0, γ, start_u, L, use_reflective)
    val = 0.0
    for n in 100000:-1:0
        if use_reflective
            val += exp(-((n * pi / L)^2) * γ / 4) * reflective_bound_exact(x, L, n) * reflective_bound_exact(x0, L, n)
        else
            val += exp(-((n * pi / L)^2) * γ / 4) * absorbing_bound_exact(x, L, n) * absorbing_bound_exact(x0, L, n)
        end
    end
    return val * start_u
end

function reflective_bound_exact(x, L, n)
    if n == 0
        return sqrt(1 / L)
    else
        return sqrt(2 / L) * cos(n * pi * x / L)
    end
end

function absorbing_bound_exact(x, L, n)
    if n == 0
        return 0.0
    else
        return sqrt(2 / L) * sin(n * pi * x / L)
    end
end

function gaussian_unbound(x, x0, γ, start_u)
    return (start_u / sqrt(γ * pi)) * exp(-(x - x0)^2 / γ)
end

function check_conservation(u, N)
    println(sum(u[1:N]) / (N + 1))
end

function eval_D(x, x0, D_plus, D_minus)
    if x >= 0.5
        return D_plus
    else
        return D_minus
    end
end

function A_plus(t, x0, D_plus, D_minus)
    result = (D_plus - D_minus) * x0^2 / (4 * D_minus * D_plus * t)
    result = exp(result) * sqrt(D_minus / D_plus) * (1 - erf(-x0 / sqrt(4 * D_minus * t)))
    result += 1 + erf(-x0 / sqrt(4 * D_plus * t))
    return 2 / result
end

function A_minus(t, x0, D_plus, D_minus)
    result = A_plus(t, x0, D_plus, D_minus) * sqrt(D_minus / D_plus)
    return result * exp((D_plus - D_minus) * x0^2 / (4 * D_minus * D_plus * t))
end

function init_u(x)
    return max(1-16*(x-0.5)^2,0)
end

function update_u(u, γ)
    updated_u = (1-γ^2)*u
    updated_u[1:N-1] += γ/2*(γ-1)*u[2:N]
    updated_u[N] += γ/2*(γ-1)*u[1] #This to have u(0)=u(1)
    updated_u[2:N] +=  γ/2*(γ+1)*u[1:N-1]
    updated_u[1] +=  γ/2*(γ+1)*u[N] #This to have u(0)=u(1)

    return updated_u
end

function init_u_gauss(x)
    σ = 0.02
    return ( sqrt(σ * pi)) * exp(-(x - 0.5)^2 / σ)
end

function oldupdate_u_hopf(u, γ)
    updated_u = u
    for n in 2:N-1
        updated_u[n] -= γ/4*(u[n+1]^2-u[n-1]^2)  
        updated_u[n] += γ^2/8*(u[n+1]+u[n])*(u[n+1]^2-u[n]^2)
        updated_u[n] -= γ^2/8*(u[n]+u[n-1])*(u[n]^2-u[n-1]^2) 
    end

    return updated_u
end

function update_u_hopf(u, γ)
    updated_u = u
    for n in 2:N-1
        updated_u[n] -= γ/4*(u[n+1]^2-u[n-1]^2)  
        updated_u[n] += γ^2/16*(u[n+1])*(u[n+1]^2-u[n]^2)
        updated_u[n] -= γ^2/16*(u[n-1])*(u[n]^2-u[n-1]^2) 
    end

    return updated_u
end


dt = 0.01
N = 85
x = [i / (N + 1) for i in 1:N]
@show dx = 1 / (N + 1)
@show γ = dt / dx

n_steps = 500

u = [init_u_gauss(x_i) for x_i in x]

anim = @animate for k in 1:n_steps
    global u = update_u_hopf(u,γ)
    if k%1==0
        plot(x, u, ylim=(0.0, 0.3), label="Lax-Wendroff", title="Step: $k", show=true)
    end
end
gif(anim, "/home/frossi/ComputationalPhysics/Assignment_1/Lax_wendroff.gif", fps=100)

# display(plot(x, u, label="Lax-Wendroff", show=true))








# dirichlet(N, α, start_u / dx, n_steps)
# neumann(N, α, start_u / dx, n_steps)
# γ = 4 * D * dt * n_steps
# unbound = [gaussian_unbound(x_i, x0, γ, start_u) for x_i in x]

# display(plot!(x, unbound, label="Unbound_exact", show=true))

# exact_reflective = [analytical_solution(x_i, x0, γ, start_u, L, true) for x_i in x]
# plot(x, exact_reflective, label="Exact_reflective", show=true)
# exact_absorbing = [analytical_solution(x_i, x0, γ, start_u, L, false) for x_i in x]
# display(plot!(x, exact_absorbing, label="Exact_absorbing", show=true))

# dirichlet_space_dependent(N, α_vec, start_u/dx, n_steps)
# neumann_space_dependent(N, α_vec, start_u/dx, n_steps)

# unbound_sp_dependent = [start_u* unbounded_space_dependent(x_i, x0, dt*n_steps, D_plus, D_minus) for x_i in x]
# display(plot!(x, unbound_sp_dependent, label="Unbound_dependent", show=true))