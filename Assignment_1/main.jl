using LinearAlgebra

function next_step(C, u)
    # Looks like it should be already ok when writing the matrix C
    return C * u
end

function dirichlet(N, α, start_u, steps)
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
        global u_n = C * u_n
        println(u_n)
    end
end

function neumann(N, α, start_u, steps)
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
        global u_n = C * u_n
        check_conservation(u_n)
        println(u_n)
    end
end

function check_conservation(u)
    println(sum(u))
end

dt = 0.00001
N = 51
D = 1
@show dx = 1 / (N + 1)
α = D * dt / (dx * dx)
start_u = 1.0

dirichlet(N, α, start_u, 5)
neumann(N, α, start_u, 5)
