using LinearAlgebra
using SparseArrays
using Printf
using Random
using Distributions
include("frequencyPrediction.jl")
include("networkModification.jl")
using DataFrames
using ArgParse


# for two-strategy games, use multiple_donationGame_K1K2 to calculate the structure coefficients 
# for multi-strategy games, use multiple_calculate_structure_coefficient to calculate the structure coefficients 


# An example
alpha = 0
numOfStrategies = 3
mutation_rate = 0.01

# generate a scale-free network of size 50 with average degree 4
num_nodes = 50
degree = 4
edge_list = generate_and_export_power_graph(num_nodes, degree, 1)

# group partation for each node with gamma = 1.0, M = 100 (\overline{M} = 2)
M = 100
gamma = 1.0
group_list, group_data = split_edges(edge_list, num_nodes, M, :weighted_absolute_deterministic, gamma)



# calculate structure coefficients lambda1, lambda2, lambda3 under DB updating
lambda1, lambda2, lambda3= multiple_calculate_structure_coefficient(edge_list, num_nodes, group_list, group_data, numOfStrategies, mutation_rate, alpha)
println("lambda1:$(lambda1)\tlambda2:$(lambda2)\tlambda3:$(lambda3)")

# calculate structure coefficients K1, K2 under DB updating
K1, K2 = multiple_donationGame_K1K2(edge_list, num_nodes, group_list, group_data, mutation_rate, alpha)
println("K1:$(K1)\tK2:$(K2)")  

#print critical benefit-to-cost ratio
println("bc_critical: $(K1/K2)")
println("bc_critical: $((lambda1+lambda2+lambda3)/(lambda1-lambda2))")

