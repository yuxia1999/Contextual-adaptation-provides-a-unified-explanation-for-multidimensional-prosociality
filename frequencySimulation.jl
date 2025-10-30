using LinearAlgebra
using Random
using Distributions 
using ArgParse
using SparseArrays
using DelimitedFiles

function calculate_offset(group_list)
    offset = [0; cumsum(group_list)[1:end-1]]
    return offset
end

function edgeList2adjacencyMatrix(edge_list, num_nodes)
    adj_matrix = spzeros(num_nodes, num_nodes)
    for edge in edge_list
        row, col, weight = edge
        adj_matrix[row+1, col+1] = weight
    end
    return adj_matrix
end

function adj2highorder_adj(adj, num_nodes, group_list, offset, M, group_data)
    W = spzeros(M, num_nodes)
    for i in 1:num_nodes
        for m in 1:group_list[i]
            for j in 1:num_nodes
                if (j-1) in group_data[i][m]
                    W[m + offset[i], j] = adj[i, j]
                end
            end
        end
    end
    return W
end


function calculate_D(edge_list, num_nodes, group_list, group_data, offset)
    D = spzeros(sum(group_list), sum(group_list))
    for i in 1:num_nodes
        for m in 1:group_list[i]
            for j in 1:num_nodes
                for n in 1:group_list[j]
                    if  (j-1) in group_data[i][m]
                        D[m+offset[i], n+offset[j]] = 1
                    end
                end
            end
        end
    end
    return D
end

function multiple_calculate_omega(W, num_nodes, offset, group_list)
    M = sum(group_list)
    omega_matrix = zeros(M,M)

    W_sign = W .!= 0

    for i in 1:num_nodes
        for m in 1:group_list[i]
            for h in 1:num_nodes
                for s in 1:group_list[h]
                    omega_matrix[m+offset[i], s + offset[h]] = W[m+offset[i], h] * W_sign[s+offset[h], i]
                end
            end
        end
    end

    return omega_matrix
    
end

function state_initialization(M, numOfStrategies)
    matrix = zeros(Int, numOfStrategies, M)

    for j in 1:M
        i = rand(1:numOfStrategies)  
        matrix[i, j] = 1
    end
    
    return matrix
end


function calculate_fitness(u, delta)
    f = exp.(u * delta)
    return f
end


function calculate_payoff(PayoffMatrix, x, omega, num_nodes, offset, group_list, M)
    u  = zeros(M)

    u1_temp = transpose(x)*PayoffMatrix*x
    u1 = u1_temp.*omega
    u1 = sum(u1, dims=2)

    for i in 1:num_nodes
        start_idx = offset[i] + 1
        end_idx = offset[i] + group_list[i]
        u_temp = sum(u1[start_idx:end_idx])
        u[start_idx:end_idx] .= u_temp
    end
    
    return u
end


function update_state(D, f, x, mutation_rate, numOfStrategies,num_nodes, group_list, offset)
    j = rand(1:num_nodes)  # select the player to be updated
    n = rand(1:group_list[j])
    index = offset[j] + n
    e = f .* D[:, index]
    e /= sum(e)
    selected_element = rand(Categorical(vec(e)))
    if rand()<=mutation_rate
        strategy = rand(1:numOfStrategies) #mutate
        x[:,index] .= 0
        x[strategy, index] = 1
    else
        x[:,index] = x[:,selected_element]
    end
    return x
end

function simulation(edge_list, num_nodes, group_list, group_data, PayoffMatrix, numOfStrategies, delta, mutation_rate, iteration_number)
    offset = calculate_offset(group_list)
    M = sum(group_list)
    adj = edgeList2adjacencyMatrix(edge_list, num_nodes)
    W = adj2highorder_adj(adj, num_nodes, group_list, offset, M, group_data)
    D = calculate_D(edge_list, num_nodes, group_list, group_data, offset)
    omega =  multiple_calculate_omega(W, num_nodes, offset, group_list)

    x = state_initialization(M, numOfStrategies)
    strategy_count = zeros((numOfStrategies, M))

    for k in 1:1000000
        u = calculate_payoff(PayoffMatrix, x, omega, num_nodes, offset, group_list,M)
        f = calculate_fitness(u,delta)
        x = update_state(D, f, x, mutation_rate, numOfStrategies, num_nodes, group_list, offset)
    end
    for k in 1:iteration_number
        u = calculate_payoff(PayoffMatrix, x, omega, num_nodes, offset, group_list,M)
        f = calculate_fitness(u,delta)
        x = update_state(D, f, x, mutation_rate, numOfStrategies, num_nodes, group_list, offset)
        strategy_count += x
    end

    results = sum(strategy_count, dims = 2)
    return results
end


function read_edge_list(file_path::String)
    edge_list = Vector{Tuple{Int, Int, Int}}()
    open(file_path, "r") do file
        for line in eachline(file)
            parts = split(chomp(line), ',')
            if length(parts) == 3
                try
                    edge = (parse(Int, parts[1]), parse(Int, parts[2]), parse(Int, parts[3]))
                    push!(edge_list, edge)
                catch e
                    println("Error parsing line: $line")
                    println(e)
                end
            else
                println("Invalid line format: $line")
            end
        end
    end
    return edge_list
end


# Example usage:
# It may takes 6 to 8 hours to finish one simulation with the following parameters.

# num_nodes = 50 
# delta = 0.005
# iteration_number = 10000000
# file_path = "SF_50_degree4_1.txt"
# edge_list = read_edge_list(file_path)
# numOfStrategies = 2
# mutation_rate = 0.01
# benefit = 5.0
# cost = 1.0    
# group_list = [5, 4, 6, 7, 5, 3, 2, 3, 2, 2, 3, 2, 2, 2, 1, 2, 1, 4, 2, 2, 1, 2, 2, 1, 2, 2, 3, 2, 2, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1]
# group_data = [[[14, 5, 39], [3, 8], [2, 42], [15, 43], [4, 1]], [[4, 0, 10], [2, 26, 28], [13, 3], [21, 39]], [[4, 27, 17], [3, 9, 46], [34, 10, 6], [7, 13], [1, 5], [16, 0]], [[44, 20, 4], [18, 6, 11], [32, 2, 22], [42, 14, 28], [29, 0], [26, 12], [1, 25]], [[36, 1, 7], [3, 2], [30, 40], [0, 12], [19, 38]], [[31, 8], [20, 2], [0]], [[2, 16], [3, 23]], [[23, 9, 4], [37, 2], [15, 11]], [[19, 5], [0]], [[7, 2], [49]], [[40, 25], [18, 2], [1, 30]], [[3, 17], [7]], [[24, 3], [4]], [[1], [2]], [[3, 0]], [[35, 7], [32, 0]], [[2, 6]], [[11, 33], [31, 2], [22, 43], [48, 36]], [[38, 3], [10]], [[24, 8], [21, 4]], [[5, 3]], [[1, 19], [35, 46]], [[33, 3], [29, 17]], [[7, 6]], [[19, 47], [12]], [[10, 27], [3]], [[34, 47], [3, 41], [48, 1]], [[25, 41], [2]], [[44, 3], [1]], [[3, 22]], [[4, 10]], [[17, 5]], [[3, 15]], [[17, 22]], [[26, 2]], [[15, 45], [21, 37]], [[17, 4]], [[35, 7]], [[18, 4]], [[0, 1]], [[10, 4]], [[26, 49], [45, 27]], [[0, 3]], [[17, 0]], [[3, 28]], [[35, 41]], [[21, 2]], [[24, 26]], [[17, 26]], [[41, 9]]]
# PayoffMatrix = [(benefit - cost) -cost; benefit 0]
# num = simulation(edge_list, num_nodes, group_list, group_data, PayoffMatrix, numOfStrategies, delta, mutation_rate, iteration_number)
# println("Frequency of strategy 1: ")
# println(num[1]/sum(num))
