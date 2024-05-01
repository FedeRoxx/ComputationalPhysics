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

function compare_in_plot(v_lists, labels, timesteps)
    for k in timesteps # 
        display(plot_list(v_lists, labels, k))
    end
end


#########################
# Problem 3, first part #
#########################

dt = 0.00015 # for α=0.375
# dt = 0.0002004  # for α=0.501
# dt = 0.001 # for α=2.5

L = 1.0  # Interval is [0,L]
N = 51 # n of points
@show dx = L/(N-1)
x = [i*dx for i in 0:N-1]
x0 = L/2
λ = 1.0
τ = 1.0
β = dt / τ
@show α = λ^2*β/dx^2

n_steps = 100

if false   # Set this to true to run 
    v = initialize_gaussian(x, x0, 0.01)
    @show V0 = trapezoidal_integral(v, dx)
    v = v / V0
    display(plot_potential(x, v, "Initial potential"))

    v_list_EE = run_euler_explicit(v, α, β, n_steps)
    v_list_EI = run_euler_implicit(v, α, β, n_steps)
    v_list_CN = run_CN(v, α, β, n_steps)
    v_list_exact = run_exact_unbound(x, x0, λ, τ, 1.0, dt, n_steps)

    compare_with_gif([v_list_EE, v_list_EI, v_list_CN, v_list_exact], ["Explicit euler", "Implicit euler", "Crank Nicolson", "Exact unbound"])
    compare_in_plot([v_list_EE, v_list_EI, v_list_CN, v_list_exact], ["Explicit euler", "Implicit euler", "Crank Nicolson", "Exact unbound"], [25,100,500])
end

##########################
# Problem 3, second part #
##########################

function g_Na(V, V_star, γ, g_K)
    g_Na = exp(γ*(V_star - V))
    g_Na = 100/(1+g_Na) + 0.2
    return g_Na/g_K
end

function CN_iteration(v, α, β, V_Na, V_K, γ, V_star, g_K)
    # Get C, q to solve v_n+1 = C v_n + q
    N = length(v)

    g_na_vec = [g_Na(V_i, V_star, γ, g_K) for V_i in v]
    diag_plus = ones(N) * (1 + α + β/2)
    diag_plus +=  β/2 * g_na_vec
    diag_minus = ones(N) * (1 - α - β/2)
    diag_minus -=  β/2 * g_na_vec

    upp_diag = -α/2* ones(N-1)
    upp_diag[1] = -α
    low_diag = -α/2* ones(N-1)
    low_diag[end] = -α
    A = Tridiagonal(low_diag, diag_plus, upp_diag)
    B = Tridiagonal(-low_diag, diag_minus, -upp_diag)

    q = ones(N) * β * V_K
    q += g_na_vec * β * V_Na
    v = B*v + q
    v = inv(A) * v

    return v
end

function EE_iteration(v, α, β, V_Na, V_K, γ, V_star, g_K)
    # Get C, q to solve v_n+1 = C v_n + q
    N = length(v)

    g_na_vec = [g_Na(V_i, V_star, γ, g_K) for V_i in v]
    diag_minus = ones(N) * (1 - 2α - β)
    diag_minus -=  β * g_na_vec

    upp_diag = -α* ones(N-1)
    upp_diag[1] = -2α
    low_diag = -α* ones(N-1)
    low_diag[end] = -2α
    B = Tridiagonal(-low_diag, diag_minus, -upp_diag)

    q = ones(N) * β * V_K
    q += g_na_vec * β * V_Na
    v = B*v + q

    return v
end

function run_CN_channels(v, α, β, n_steps, V_Na, V_K, γ, V_star, g_K)
    v_list = []
    anim = @animate for k in 1:n_steps
        # v_old = deepcopy(v)
        v = CN_iteration(v, α, β, V_Na, V_K, γ, V_star, g_K)
        # if 10<k<30
        #     @show g_na_vec_old = [g_Na(V_i, V_star, γ, g_K) for V_i in v]
        #     @show g_na_vec = [g_Na(V_i, V_star, γ, g_K) for V_i in v_old]
        # end
        plot_potential_channel(x, v, "Crank Nicolson with "*L"V_{appl}="*"-40 mV, step "*string(k))
        push!(v_list, v)
    end
    display(gif(anim, "/home/frossi/ComputationalPhysics/Exam/Prob3/Crank_Nicolson_channel.gif", fps=10))
    display(plot_potential_channel(x, v, "Crank Nicolson with "*L"V_{appl}="*"-40 mV"))
    return v_list
end

function run_EE_channels(v, α, β, n_steps, V_Na, V_K, γ, V_star, g_K)
    v_list = []
    anim = @animate for k in 1:n_steps
        v = EE_iteration(v, α, β, V_Na, V_K, γ, V_star, g_K)
        plot_potential_channel(x, v, "Euler explicit, step "*string(k))
        push!(v_list, v)
    end
    display(gif(anim, "/home/frossi/ComputationalPhysics/Exam/Prob3/Euler_explicit_channel.gif", fps=20))
    display(plot_potential_channel(x, v, "Euler explicit"))
    return v_list
end

function initialize_gauss_channel(x, V_appl, V_mem, x0, λ)
    gaussian = exp.(-0.5 * ((x .- x0) / λ).^2) 
    return (V_appl-V_mem)*gaussian .+ V_mem
end

function plot_potential_channel(x, v, title)
    plot(x, v, ylim=(-75, 50), xlabel=L"x"*" [mm]", ylabel=L"V(x)"*" [mV]", title=title, legend=false)
end

function plot_list_channel(v_lists, labels, k)
    p = plot(ylim=(-75, 50), label=labels[1], xlabel=L"x", ylabel=L"V(x)", title="All methods compared, step "*string(k), legend=true)
    for (i, v_list) in enumerate(v_lists)
        if labels[i] == "Exact unbound"
            plot!(p, x, v_list[k], label=labels[i])
        else
            plot!(p, x, v_list[k], label=labels[i], marker=:x)
        end
    end
    return plot!(p)
end

function compare_with_gif_channel(v_lists, labels)
    n_steps = length(v_lists[1])
    anim = @animate for k in 1:n_steps # 
        plot_list_channel(v_lists, labels, k)
    end
    display(gif(anim, "/home/frossi/ComputationalPhysics/Exam/Prob3/Comparison_channel.gif", fps=10))
end

function plot_potential_with_time(v_list, point, dt)
    N = length(v_list)
    t = [dt*i for i in 1:N]
    v = [v_i[point] for v_i in v_list]
    plot(t, v, ylim=(-75, 50), xlabel=L"t", ylabel=L"V(t)", title="Potential at some position", legend=false)
end

dt = 0.002 # for α=0.375
# dt = 0.0002004  # for α=0.501
# dt = 0.001 # for α=2.5

L = 1.0  # Interval is [0,L]    
N = 51 # n of points
@show dx = L/(N-1)
x = [i*dx for i in 0:N-1]
x0 = L/2
λ = 0.18
τ = 2.0
β = dt / τ
@show α = λ^2*β/dx^2
γ = 0.5
V_star = -40
V_Na = 56
V_K = -76
g_K = 5.0

###NOTE###
# Now times are in ms, lengths in mm and V in mV

n_steps = 500
V_appl = -50
V_mem = -70

### Run with V_appl = -50 [Task a]
v = initialize_gauss_channel(x, V_appl, V_mem, x0, λ)
v_list_channel = run_CN_channels(v, α, β, n_steps, V_Na, V_K, γ, V_star, g_K)
display(plot_potential_with_time(v_list_channel, 30, dt))

### Testing different V_appl [Task b]
p = plot(ylim=(-75, 50), xlabel=L"t"*" [ms]", ylabel=L"V(t)"*" [mV]", title="Evolution of potential at "*L"x=0.76"*" mm", legend=:right)
for V_appl in [-60, -50, -47, -46, -43, -40]
    v = initialize_gauss_channel(x, V_appl, V_mem, x0, λ)
    new_v_list_channel = run_CN_channels(v, α, β, n_steps, V_Na, V_K, γ, V_star, g_K)
    t = [dt*i for i in 1:n_steps]
    v_points = [v_i[38] for v_i in new_v_list_channel]
    plot!(p, t, v_points, label=L"V_{appl}="*string(V_appl)*" mV")
end
display(plot!(p))


# ###Testing Euler explicit
# v = initialize_gauss_channel(x, V_appl, V_mem, x0, λ)
# v_list_channel_EE = run_EE_channels(v, α, β, n_steps, V_Na, V_K, γ, V_star, g_K)
# compare_with_gif_channel([v_list_channel, v_list_channel_EE], ["CN", "EE"])


##############################
# Problem 3, shifted channel #
##############################

function CN_iteration_shifted(v, α, β, V_Na, V_K, γ, V_star, g_K)
    # Get C, q to solve v_n+1 = C v_n + q
    N = length(v)

    # This time g is kept zero up to 1mm
    g_na_vec = [g_Na(V_i, V_star, γ, g_K) for V_i in v]
    g_na_vec[1:20] .= 0.0
    # g_na_vec = zeros(N)
    # g_na_vec[21] = g_Na(v[21], V_star, γ, g_K)
    diag_plus = ones(N) * (1 + α + β/2)
    diag_plus +=  β/2 * g_na_vec
    diag_minus = ones(N) * (1 - α - β/2)
    diag_minus -=  β/2 * g_na_vec

    upp_diag = -α/2* ones(N-1)
    upp_diag[1] = -α
    low_diag = -α/2* ones(N-1)
    low_diag[end] = -α
    A = Tridiagonal(low_diag, diag_plus, upp_diag)
    B = Tridiagonal(-low_diag, diag_minus, -upp_diag)

    q = ones(N) * β * V_K
    q += g_na_vec * β * V_Na
    v = B*v + q
    v = inv(A) * v

    return v
end

function run_CN_channels_shifted(v, α, β, n_steps, V_Na, V_K, γ, V_star, g_K)
    v_list = []
    anim = @animate for k in 1:n_steps
        v_old = deepcopy(v)
        v = CN_iteration_shifted(v, α, β, V_Na, V_K, γ, V_star, g_K)
        if 10<k<30
            @show g_na_vec_old = [g_Na(V_i, V_star, γ, g_K) for V_i in v]
            @show g_na_vec = [g_Na(V_i, V_star, γ, g_K) for V_i in v_old]
        end
        plot_potential_channel(x, v, "Crank Nicolson, step "*string(k))
        push!(v_list, v)
    end
    display(gif(anim, "/home/frossi/ComputationalPhysics/Exam/Prob3/Crank_Nicolson_shifted_channel.gif", fps=10))
    display(plot_potential_channel(x, v, "Crank Nicolson"))
    return v_list
end

dt = 0.01 # for α=0.375
# dt = 0.0002004  # for α=0.501
# dt = 0.001 # for α=2.5

L = 2.5  # Interval is [0,L]    
N = 51 # n of points
@show dx = L/(N-1)
x = [i*dx for i in 0:N-1]
x0 = 0.75
λ = 0.18
τ = 2.0
β = dt / τ
@show α = λ^2*β/dx^2
γ = 0.5
V_star = -40
V_Na = 56
V_K = -76
g_K = 5.0

###NOTE###
# Now times are in ms, lengths in mm and V in mV

n_steps = 500
V_appl = 0
V_mem = -70
v = initialize_gauss_channel(x, V_appl, V_mem, x0, λ)
# v_list_channel = run_CN_channels(v, α, β, n_steps, V_Na, V_K, γ, V_star, g_K)