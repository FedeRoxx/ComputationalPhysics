using LinearAlgebra
using Plots
using SpecialFunctions

function initialize_wave(N)
    #Not sure if it should be N
    u_0 = [sin(pi * i / (N)) * sin(2 * pi * j / (N)) for (i, j) in zip(0:N, 0:N)]
    return reshape(u_0, N+1, N+1)
end

function set_boundaries(u)
    u[1,:] = 0.0
    u[:,1] = 0.0
    u[end,:] = 0.0
    u[:,end] = 0.0
    return u
end

function update_u(u, u_old, β, N)
    new_u = (1-2*β)*2*u - u_old
    for (x, y) in zip(2:N, 2:N)
        new_u[x, y] += β *( u[x+1,y] + u[x-1,y] + u[x,y+1] + u[x,y-1])
    end
    return set_boundaries(new_u)
end


function exact_solution(N, t)
    u = [cos(sqrt(5)*pi*t)*sin(pi * i / (N)) * sin(2 * pi * j / (N)) for (i, j) in zip(0:N, 0:N)]
    return reshape(u, N+1, N+1)
end

dt = 0.00001
N = 51
c = 1
x = [i / N for i in 0:N]
@show dx = 1 / N
β = c^2 * dt^2 / dx^2
u = initialize_wave(N)
u = set_boundaries(u)
old_u = u
n_steps = 500


for k in 1:n_steps
    new_u = update_u(u, u_old, β, N)
    u_old = u
    u = new_u
end

u_exact = exact_solution(N, dt*n_steps)

# dirichlet(N, α, start_u / dx, n_steps)
# neumann(N, α, start_u / dx, n_steps)
# γ = 4 * D * dt * n_steps
# unbound = [gaussian_unbound(x_i, x0, γ, start_u) for x_i in x]

# display(plot!(x, unbound, label="Unbound_exact", show=true))

# exact_reflective = [analytical_solution(x_i, x0, γ, start_u, L, true) for x_i in x]
# display(plot(x, exact_reflective, label="Exact_reflective", show=true))
# exact_absorbing = [analytical_solution(x_i, x0, γ, start_u, L, false) for x_i in x]
# display(plot(x, exact_absorbing, label="Exact_absorbing", show=true))

# dirichlet_space_dependent(N, α_vec, start_u, n_steps)
# neumann_space_dependent(N, α_vec, start_u, n_steps)

# unbound_sp_dependent = [start_u / dx * unbounded_space_dependent(x_i, x0, dt * n_steps, D_plus, D_minus) for x_i in x]
# display(plot!(x, unbound_sp_dependent, label="Unbound_dependent", show=true))
