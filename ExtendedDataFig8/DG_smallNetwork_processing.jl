using Printf
using Random
using Distributions
include("frequencyPrediction.jl")
include("networkModification.jl")
using DataFrames
using ArgParse
using DelimitedFiles


function DG_bcCritical(lambda1, lambda2, lambda3, selection_intensity)

    return (lambda1+lambda2+lambda3)/(lambda1-lambda2)
end



selection_intensity = 0.001
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

bc_ori = DG_bcCritical(lambda1_ori, lambda2_ori, lambda3_ori, selection_intensity)
bc_grouped = DG_bcCritical(lambda1_grouped, lambda2_grouped, lambda3_grouped, selection_intensity)
open("donationGame_networkOfSize6_bc_ungrouped_u0.1.txt", "a") do io
        println(io, @sprintf("%.10f", bc_ori))
end

open("donationGame_networkOfSize6_bc_grouped_u0.1.txt", "a") do io
    println(io, @sprintf("%.10f", bc_grouped))
end

end
