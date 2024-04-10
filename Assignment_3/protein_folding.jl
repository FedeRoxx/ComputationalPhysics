using LinearAlgebra
using Plots
using SpecialFunctions
using LaTeXStrings
using Printf
# using Optim
using Serialization
using Random

global dict_direction = Dict(
    1 => [1, 0],
    2 => [-1, 0],
    3 => [0, 1],
    4 => [0, -1]
)

function add_2d_monomer(chain)
    # Adds a monomer at the end of the protein, with random type and dict_direction
    # If no positions are available, builds a new chain of same length from scratch
    if length(chain)==1
        push!(chain, [rand(1:20), dict_direction[rand(1:4)]])
    else
        shuffled_list = shuffle([1,2,3,4])
        found_dir = false
        for i in shuffled_list
            new_dir = dict_direction[i]
            if chain[end][2] + new_dir ∉ [point[2] for point in chain]
                push!(chain, [rand(1:20), chain[end][2] + new_dir])
                found_dir = true
                break
            end
        end
        if !found_dir
            println("Random generation stuck! Restarted from scratch")
            chain = build_2d_protein(length(chain)+1)
        end
    end
    return chain
end

function build_2d_protein(n)
    # Builds a protein with n monomers of random types, starting in 0,0
    # Protein is respresented as a list of monomers [type_i, [x_i, y_i]] with type_i in 1:20
    protein = [[rand(1:20), [0, 0]]]
    for _ in 1:n-1
        protein = add_2d_monomer(protein)
    end
    return protein
end

function plot_2d_protein(protein, title)
    # Plots a 2d graph with monomers indicated by their type and connected with blue lines
    x_coords = [point[2][1] for point in protein]
    y_coords = [point[2][2] for point in protein]
    labels = [point[1] for point in protein]

    # scatter(x_coords, y_coords, xlims=(-5,20), ylims=(-10,10), label = "")
    scatter(x_coords, y_coords, aspect_ratio=:equal, label = "")

    # for i in 1:length(protein)
    #     annotate!(x_coords[i], y_coords[i], text(labels[i], 15, :black))
    # end
    for i in 1:length(protein) - 1
        plot!([x_coords[i], x_coords[i+1]], [y_coords[i], y_coords[i+1]], color=:blue, linewidth=2, label="")
    end
    xlabel!("X")
    ylabel!("Y")
    title!(title)

    display(plot!())
end

function get_J_mat(filename)
    # If already present, reads the J_mat saved.
    # Otherwise, builds a 20x20 symmetric matrix of random values between -2 and -4
    if isfile(filename)
        J_mat = deserialize(filename)
    else
        J_mat = -2 .- 2 .* rand(20, 20)
        for i in 1:20, j in i+1:20
            J_mat[j, i] = J_mat[i, j]
        end
        serialize(filename, J_mat)
    end
    return J_mat
end

function generate_nearest_neighbour(protein)
    # For each monomer, builds a list of the types of the nearest neighbours
    NN_list = []
    for n1 in eachindex(protein)
        local_NN = []
        for n2 in eachindex(protein)
            if any(x -> x == protein[n2][2] - protein[n1][2], values(dict_direction)) && abs(n1-n2)>1
                push!(local_NN, protein[n2][1])
            end
        end
        push!(NN_list, local_NN)
    end
    return NN_list           
end

function evaluate_energy(protein, J_mat)
    # Builds the neighbours list and evaluates energy, with a given J_mat
    NN_list = generate_nearest_neighbour(protein)
    energy = 0.0
    for (i, NN) in enumerate(NN_list)
        for el in NN
            energy += J_mat[protein[i][1], el]
        end
    end
    return energy
end

function find_move(index, protein)
    # Find a random available move for the monomers in position index.
    # Returns true, new_coord (=the new suggested coordinates for the monomer)
    # If none are available returns false, [0,0]
    coords = [monomer[2] for monomer in protein]
    shuffled_list = shuffle([1,2,3,4])
    if index == 1
        for i in shuffled_list
            new_coord = coords[index+1] + dict_direction[i]
            if new_coord ∉ coords
                return true, new_coord
            end
        end
    elseif index == length(protein)
            for i in shuffled_list
                new_coord = coords[index-1] + dict_direction[i]
                if new_coord ∉ coords
                    return true, new_coord
                end
            end
    else
        new_coord = coords[index-1]+coords[index+1]-coords[index]
        if new_coord ∉ coords
            return true, new_coord
        end
    end
    return false, [0,0]
end

function MC_step(protein, J_mat, T)
    # Complete a MC step looking for a random move on a random monomer.
    # Accept of reject the transition according to Metropolis
    energy1 = evaluate_energy(protein, J_mat)
    protein2 = deepcopy(protein)
    valid_step = false
    shuffled_list = shuffle(1:length(protein))
    for i in shuffled_list
        bool, new_coord = find_move(i, protein)
        if bool
            protein2[i][2] = new_coord
            valid_step = true
            break
        end
    end
    if !valid_step
        #This should never happen
        println("No moves are available, impossible!")
    end
    energy2 = evaluate_energy(protein2, J_mat)
    P_accept = min(1.0, exp(-(energy2-energy1)/T))
    if rand() < P_accept
        return protein2
    end
    return protein
end

function MC_sweep(protein, J_mat, T)
    # Performs N MonteCarlo steps, completing a sweep
    for _ in 1:length(protein)
        protein = MC_step(protein, J_mat, T)
    end
    @show energy = evaluate_energy(protein, J_mat)
    return protein
end

function MC_simulation(protein, J_mat, T, n_sweeps)
    # Run a simulations of n_sweeps from a starting protein.
    # Keeps a logger to get the energy, e2e and RoG.
    logger = []
    for k in 1:n_sweeps
        protein = MC_sweep(protein, J_mat, T)
        if k in [1,10,100]
            plot_2d_protein(protein, "2D representation of protein after $k sweeps")
        end
        push!(logger, get_info(protein, J_mat))
    end
    return logger
end

function CoM(protein)
    # Evaluates the Center of mass of a protein, assuming m=1 for each monomer.
    coord = [0,0]
    for monomer in protein
        coord += monomer[2]
    end
    return coord / length(protein)
end

function RoG(protein)
    # Calculate Radius of Gyration for a given protein
    com = CoM(protein)
    temp_sum = sum([norm(monomer[2] - com)^2 for monomer in protein])
    return sqrt(temp_sum/length(protein))
end

function get_info(protein, J_mat)
    # For a gived protein, evaluates energy, e2e and RoG and give them back.
    info = [evaluate_energy(protein, J_mat)]
    push!(info, norm(protein[end][2]-protein[1][2]))
    push!(info, RoG(protein))
    return info
end

function save_protein(protein, filename)
    # Saves the given protein in the filename.
    serialize(filename, protein)
    return nothing
end

function load_protein(filename)
    # Loads a protein from the file given as filename.
    return deserialize(filename)
end

function linear_protein(n)
    protein = [[rand(1:20), [i, 0]] for i in 0:n-1]
    return protein
end

function display_logger(logger)
    # Plots energy, e2e and RoG as a function of the sweeps index.
    energies = [point[1] for point in logger]
    e2e = [point[2] for point in logger]
    RoG = [point[3] for point in logger]

    display(plot(1:length(logger), energies, label="Energy", xlabel="Sweeps", linewidth=1, title="Energy during the simulation"))
    display(plot(1:length(logger), e2e, label="End-to-end", xlabel="Sweeps", linewidth=1, title="End-to-end during the simulation"))
    display(plot(1:length(logger), RoG, label="RoG", xlabel="Sweeps", linewidth=1, title="Radius of gyration during the simulation"))
    return nothing
end

function MC_simulation_averaged(protein, J_mat, T, n_sweeps, n_average)
    # Run a simulations of n_sweeps from a starting protein.
    # Then returns energy, e2e and RoG averaged for more n_average sweeps.
    logger = [[0.0,0.0,0.0]]
    plot_2d_protein(protein, "Before starting")
    for k in 1:n_sweeps
        protein = MC_sweep(protein, J_mat, T)
        push!(logger, get_info(protein, J_mat))
        if k % 200 == 0
            plot_2d_protein(protein, "2D protein at "*L"T=1")
        end
    end
    means = [0.0,0.0,0.0]
    for k in 1:n_average
        protein = MC_sweep(protein, J_mat, T)
        means += get_info(protein, J_mat)
        push!(logger, means/k)
        # means += get_info(protein, J_mat)
    end
    plot_2d_protein(protein, "2D representation of protein")
    display_logger(logger)
    return means / n_average, protein
end


function MC_simulation_Tlist(protein, J_mat, T_list, n_sweeps)
    # Run a simulations of n_sweeps from a starting protein, changing to the next T in T_list after.
    # Keeps a logger to get the energy, e2e and RoG.
    logger = []
    for temp in T_list
        for _ in 1:n_sweeps
            protein = MC_sweep(protein, J_mat, temp)
            push!(logger, get_info(protein, J_mat))
        end
    end
    return logger
end

function MC_simulation_Tlist_averaged(protein, J_mat, T_list, n_sweeps, n_average)
    # Run a simulations of n_sweeps from a starting protein, changing to the next T in T_list after.
    # Keeps a logger to get the energy, e2e and RoG.
    logger = []
    for temp in T_list
        mean, protein = MC_simulation_averaged(protein, J_mat, temp, n_sweeps, n_average)
        push!(logger, mean)
    end
    return logger
end

function display_logger_Tlist(logger, T_list)
    # Plots energy, e2e and RoG as a function of the sweeps index.
    energies = [point[1] for point in logger]
    e2e = [point[2] for point in logger]
    RoG = [point[3] for point in logger]

    # display(plot(T_list, energies, label="Energy", xlabel="Temperature", linewidth=1, title="Energy at different temperatures"))
    # display(plot(T_list, e2e, label="End-to-end", xlabel="Temperature", linewidth=1, title="End-to-end at different temperatures"))
    # display(plot(T_list, RoG, label="RoG", xlabel="Temperature", linewidth=1, title="Radius of gyration at different temperatures"))
    savefig(plot(T_list, energies, label="Energy", xlabel="Temperature", linewidth=1, title="Energy at different temperatures"), "/home/frossi/ComputationalPhysics/Assignment_3/energy.png")
    savefig(plot(T_list, e2e, label="End-to-end", xlabel="Temperature", linewidth=1, title="End-to-end at different temperatures"), "/home/frossi/ComputationalPhysics/Assignment_3/e2e.png")
    savefig(plot(T_list, RoG, label="RoG", xlabel="Temperature", linewidth=1, title="Radius of gyration at different temperatures"), "/home/frossi/ComputationalPhysics/Assignment_3/rog.png")
    return nothing
end

function MC_simulation_SA(protein, J_mat, T1, T2, n_sweeps)
    # Run a simulations of n_sweeps from a starting protein, changing to the next T in T_list after.
    # Keeps a logger to get the energy, e2e and RoG.
    logger = [[0.0,0.0,0.0]]
    T_list = range(T1, stop=T2, length=n_sweeps)
    for temp in T_list
        protein = MC_sweep(protein, J_mat, temp)
        push!(logger, get_info(protein, J_mat))
    end
    plot_2d_protein(protein, "2D representation of protein")
    return logger
end



J_mat = get_J_mat("/home/frossi/ComputationalPhysics/Assignment_3/J_mat_serialized")
display(J_mat)


n = 15

if false
    protein = build_2d_protein(n)
    plot_2d_protein(protein, "2D representation of protein")
    @show NN_list = generate_nearest_neighbour(protein)
    @show energy = evaluate_energy(protein, J_mat)
    save_protein(protein, "/home/frossi/ComputationalPhysics/Assignment_3/protein_serialized" )
else
    protein = load_protein("/home/frossi/ComputationalPhysics/Assignment_3/protein_serialized")
    # plot_2d_protein(protein, "2D representation of protein")
end

### Generate and save a new linear protein of fixed length
# protein = linear_protein(30)
# save_protein(protein, "/home/frossi/ComputationalPhysics/Assignment_3/30_linear_protein_serialized" )


# 100 sweeps 
# protein = load_protein("/home/frossi/ComputationalPhysics/Assignment_3/15_linear_protein_serialized" )
# plot_2d_protein(protein, "2D representation of protein")
# @show logger = MC_simulation(protein, J_mat, 10, 100)
# display_logger(logger)

### Averaging after 100 sweeps, for 1000 sweeps [Task 5]
# @show mean, protein = MC_simulation_averaged(protein, J_mat, 10, 100, 1000)

### Temperature of 1 [Task 6 (part 1)]
# plot_2d_protein(protein, "Initial geometry of protein")
# logger = MC_simulation(protein, J_mat, 1, 100)
# display_logger(logger)

### N = 100 [Task 6 (part 2)]
# protein = load_protein("/home/frossi/ComputationalPhysics/Assignment_3/100_linear_protein_serialized" )
# plot_2d_protein(protein, "Initial geometry of protein")
# logger, protein = MC_simulation_averaged(protein, J_mat, 10, 100, 1000)
# display_logger(logger)

### Decrease temperature [Task 7 a]
# protein = load_protein("/home/frossi/ComputationalPhysics/Assignment_3/15_linear_protein_serialized" )
# T_list = [20, 10, 5, 3, 2, 1]
# plot_2d_protein(protein, "Initial geometry of protein")
# logger = MC_simulation_Tlist(protein, J_mat, T_list, 1000)
# display_logger(logger)


### Decrease temperature and take average [Task 7 b]
# protein = load_protein("/home/frossi/ComputationalPhysics/Assignment_3/50_linear_protein_serialized" )
# T_list = [20, 15, 10, 8, 6, 5, 3, 2, 1]
# plot_2d_protein(protein, "Initial geometry of protein")
# logger = MC_simulation_Tlist_averaged(protein, J_mat, T_list, 100, 3000)
# display_logger_Tlist(logger, T_list)

### Get 2 snapshots [Task 8]
# protein = load_protein("/home/frossi/ComputationalPhysics/Assignment_3/30_linear_protein_serialized" )
# mean, protein = MC_simulation_averaged(protein, J_mat, 1, 10000, 1)
# logger = MC_simulation_SA(protein, J_mat, 10.0, 1.0, 10000)
# display_logger(logger)


display(plot_2d_protein(protein, "2D representation of protein"))
