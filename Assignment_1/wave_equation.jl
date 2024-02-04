using LinearAlgebra
using Plots
using SpecialFunctions

function initialize_wave(N)
    #Not sure if it should be N
    u_x = [sin(pi * i / (N)) for i in 0:N]
    u_y = [sin(2 * pi * i / (N)) for i in 0:N]
    return u_x .* u_y'
end

function set_boundaries(u)
    u[1, :] .= 0.0
    u[:, 1] .= 0.0
    u[end, :] .= 0.0
    u[:, end] .= 0.0
    return u
end

function update_u(u, u_old, β, N)
    new_u = (1 - 2 * β) * 2 * u - u_old
    for x in 2:N
        for y in 2:N
            new_u[x, y] += β * (u[x+1, y] + u[x-1, y] + u[x, y+1] + u[x, y-1])
        end
    end
    return set_boundaries(new_u)
end


function exact_solution(N, t)
    return cos(sqrt(5) * pi * t) * initialize_wave(N)
end

dt = 0.001
N = 51
c = 1
x = [i / N for i in 0:N]
@show dx = 1 / N
β = c^2 * dt^2 / dx^2
u = initialize_wave(N)
u = set_boundaries(u)
u_old = u
n_steps = 1000


anim = @animate for k in 1:n_steps
    global u, u_old
    new_u = update_u(u, u_old, β, N)
    u_old = u
    u = new_u
    u_exact = exact_solution(N, dt * k)
    surface(u, c=:blues, zlim=(-1, 1), cbar=false, xlabel="X", ylabel="Y", title="Time Step: $k")
end
gif(anim, "/home/frossi/ComputationalPhysics/Assignment_1/wave_exact.gif", fps=100)

# display(u)

# u_exact = exact_solution(N, dt*n_steps)
# display(u_exact)
