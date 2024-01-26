using LinearAlgebra
using Plots

function next_step(C, u)
    # Looks like it should be already ok when writing the matrix C
    return C * u
end

function dirichlet(N, α, start_u, steps)
    x = [i/(N+1) for i in 1:N]
    u_n = zeros(N)
    u_n[Int((N + 1) / 2)] = start_u
    println(u_n)
 
    alpha_halved = ones(N-1) * α / 2
    alpha_plus = ones(N) * (1 + α)
    alpha_minus = ones(N) * (1 - α)
    
    A = SymTridiagonal(alpha_plus, -alpha_halved)
    B = SymTridiagonal(alpha_minus, alpha_halved)
    C = inv(A) * B
    
    for _ in 1:steps
        u_n = C * u_n
        check_conservation(u_n,N)
        println(u_n)
    end
    plot(x,u_n, label="Dirichlet", show=true)
end

function neumann(N, α, start_u, steps)
    x = [i/(N+1) for i in 1:N]
    u_n = zeros(N)
    u_n[Int((N + 1) / 2)] = start_u
    println(u_n)
 
    alpha_plus = ones(N) * (1 + α)
    alpha_minus = ones(N) * (1 - α)
    alpha_halved_1 = ones(N-1) * α / 2
    alpha_halved_1[1] = α
    alpha_halved_2 = ones(N-1) * α / 2
    alpha_halved_2[end] = α
    
    A = Tridiagonal(-alpha_halved_2, alpha_plus, -alpha_halved_1)
    B = Tridiagonal(alpha_halved_2, alpha_minus, alpha_halved_1)
    C = inv(A) * B
    
    for _ in 1:steps
        u_n = C * u_n
        check_conservation(u_n,N)
        println(u_n)
        # plot!(x,u_n, label="my label", show=true)
    end
    plot!(x,u_n, label="Neumann", show=true)
end

function exact_solution(γ, x0, start_u)
    
end

function gaussian_exact(x, x0, γ, start_u)
    return (start_u / sqrt(γ * pi)) * exp(-(x-x0)^2 / γ)
end

function check_conservation(u, N)
    println(sum(u)/(N+1))
end

dt = 0.00001
N = 51
D = 1
@show dx = 1 / (N + 1)
α = D * dt / (dx * dx)
start_u = 1.0
x0 = 0.5
n_steps = 5000

dirichlet(N, α, start_u/dx, n_steps)
neumann(N, α, start_u/dx, n_steps)
γ = 4*D*dt*n_steps
x = [i/(N+1) for i in 1:N]
exact = [gaussian_exact(x_i, x0, γ, start_u) for x_i in x]

display(plot!(x,exact, label="Exact", show=true))