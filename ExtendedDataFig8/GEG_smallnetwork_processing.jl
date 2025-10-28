using Printf
using Random
using Distributions
include("frequencyPrediction.jl")
include("networkModification.jl")
using DataFrames
using ArgParse
using DelimitedFiles

function GEG_pq_frequency_continual(lambda1, lambda2, lambda3, selection_intensity, b)
    phi = 1/2 + selection_intensity * (- lambda2 * 1/6 - lambda3 * 1/12)
    omega = 1/2 + selection_intensity * ((b-1)/12 *lambda1 - lambda2 * (b+1)/12 - lambda3 * 1/12 )
    return phi,omega
end



selection_intensity = 0.001
b = 1.5

for index in 1:112
file_path_ori = "EDF8_data/small_network_size6_index$(index)_ungrouped_structureCoefficient.txt"
file_path_grouped = "EDF8_data/small_network_size6_index$(index)_completely_grouped_structureCoefficient.txt"

data_matrix_ori = readdlm(file_path_ori, ' ')
data_matrix_grouped = readdlm(file_path_grouped, ' ')

lambda1_ori = data_matrix_ori[34,2]
lambda2_ori = data_matrix_ori[34,3]
lambda3_ori = data_matrix_ori[34,4]

lambda1_grouped = data_matrix_grouped[34,2]
lambda2_grouped = data_matrix_grouped[34,3]
lambda3_grouped = data_matrix_grouped[34,4]

p_ori,q_ori = GEG_pq_frequency_continual(lambda1_ori, lambda2_ori, lambda3_ori, selection_intensity,b)
p_grouped,q_grouped = GEG_pq_frequency_continual(lambda1_grouped, lambda2_grouped, lambda3_grouped, selection_intensity,b)
open("GEG_networkOfSize6_p_ungrouped_u0.1.txt", "a") do io
        println(io, @sprintf("%.10f", p_ori))
end
open("GEG_networkOfSize6_q_ungrouped_u0.1.txt", "a") do io
    println(io, @sprintf("%.10f", q_ori))
end
open("GEG_networkOfSize6_p_grouped_u0.1.txt", "a") do io
    println(io, @sprintf("%.10f", p_grouped))
end
open("GEG_networkOfSize6_q_grouped_u0.1.txt", "a") do io
    println(io, @sprintf("%.10f", q_grouped))
end
end
