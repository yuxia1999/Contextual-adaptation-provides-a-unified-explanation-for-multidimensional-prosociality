using LinearAlgebra
using SparseArrays
using Printf
using Random
using Distributions
using IterativeSolvers
    
function calculate_p(edge_list, num_nodes)
    edge_list = Array(edge_list)
    P = zeros(num_nodes, num_nodes)
    in_degrees = zeros(num_nodes)
    
    for edge in edge_list
        source, target, weight = edge
        in_degrees[target + 1] += weight
    end
    
    if any(in_degrees .== 0)
        println("Error: in-degree is zero")
        return nothing
    end
    
    for edge in edge_list
        source, target, weight = edge
        P[target + 1, source + 1] = weight / in_degrees[target + 1]
    end
    
    return P
end

function edgeList2adjacencyMatrix(edge_list, num_nodes)
    adj_matrix = zeros(num_nodes, num_nodes)
    
    for edge in edge_list
        row, col, weight = edge[1] + 1, edge[2] + 1, edge[3]  
        adj_matrix[row, col] = weight
    end
    
    return adj_matrix
end



function calculate_offset(group_list)
    offset = [0; cumsum(group_list)[1:end-1]]
    return offset
end

function calculate_D(edge_list, num_nodes, group_list, group_data, offset)
    D = zeros(sum(group_list), sum(group_list))
    for i in 1:num_nodes
        for m in 1:group_list[i]
            for j in 1:num_nodes
                for n in 1:group_list[j]
                    if ( j-1 in group_data[i][m])
                        D[m+offset[i], n+offset[j]] = 1
                    end
                end
            end
        end
    end
    return D
end

function calculate_A(D)
    D_transposed = transpose(D)
    row_sums = sum(D_transposed, dims=2)  
    normalized_matrix = D_transposed ./ row_sums
    
    return normalized_matrix
end

function adj2highorder_adj(adj, num_nodes, group_list, offset, group_data)
    group_list_np = Array(group_list)
    W = zeros(sum(group_list_np), num_nodes)

    for i in 1:num_nodes
        for m in 1:group_list[i]
            for j in 1:num_nodes
                if j - 1 in group_data[i][m]
                    W[m + offset[i] , j] = adj[i, j]
                end
            end
        end
    end

    return W
end

function multiple_payoffProcessing(PayoffMatrix, numOfStrategies)
    # return payoff terms

    if (size(PayoffMatrix,1) != size(PayoffMatrix,2)) || (size(PayoffMatrix,1) != numOfStrategies) || (size(PayoffMatrix,2) != numOfStrategies) 
        println("Wrong payoff matrix")
        return 0
    end

    a_yStar_bar = sum(PayoffMatrix, dims = 2) # return a column vector
    a_starY_bar = sum(PayoffMatrix, dims = 1) # return a row vector
    a_starStar_bar = 0
    a_bar = 0

    for z in 1:numOfStrategies
        a_starStar_bar += PayoffMatrix[z,z]
        for w in 1:numOfStrategies
            a_bar += PayoffMatrix[z,w]
        end
    end

    a_yStar_bar /= numOfStrategies
    a_starY_bar /= numOfStrategies
    a_starStar_bar /= numOfStrategies
    a_bar /= numOfStrategies^2

    return a_yStar_bar, a_starY_bar, a_starStar_bar, a_bar
end

function multiple_calculate_m_coefficient(A, num_nodes, offset, group_list, alpha, k,l,i,m,j,n)
    # return m_{k,l}^{(i,m),(j,n)}

    I = Diagonal(ones(num_nodes))
    sum_A_jnks = 0
    for s in 1:group_list[k]
        sum_A_jnks += A[n+offset[j], s + offset[k]]
    end

    m_kl_imjn = 1 / (num_nodes * group_list[j]) * A[n+offset[j], m+offset[i]] * ( I[i,k]*(I[l,m]*alpha + 1- alpha) - (alpha*A[n+offset[j], l+offset[k]] + (1-alpha) * sum_A_jnks) )

    return m_kl_imjn
end

function multiple_calculate_pi_mut(A, num_nodes, offset, group_list, mutation_rate)
    M = sum(group_list)
    coeff_matrix = zeros(M,M)
    b_vector = zeros(M,1)

    for i in 1:num_nodes
        for m in 1:group_list[i]
            b_vector[m+offset[i]] = -1 * mutation_rate * group_list[i] * num_nodes / M 
            coeff_matrix[m+offset[i], m+offset[i]] = -1
            for j in 1:num_nodes
                for n in 1:group_list[j]
                    coeff_matrix[m+offset[i], n+offset[j]] += (1-mutation_rate) * A[n+offset[j], m+offset[i]] * group_list[i] / group_list[j]
                end
            end
        end
    end
   
    #pi_mut_vector = coeff_matrix \ b_vector
    pi_mut_vector = idrs(coeff_matrix,b_vector) 
    
    return pi_mut_vector
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


function multiple_calculate_phi_secondOrder(A, num_nodes, offset, group_list, mutation_rate, numOfStrategies)
    M = sum(group_list)

    L_vector = Vector{Int}(undef, M)
    current_pos = 1
    for group_size in group_list
        range_end = current_pos + group_size - 1
        L_vector[current_pos:range_end] .= group_size
        current_pos = range_end + 1
    end

    nnz_estimate = round(Int, M * (M - 1) / 2 * (3 + 2 * M) + M)
    I = Vector{Int}(undef, nnz_estimate)
    J = Vector{Int}(undef, nnz_estimate)
    V = Vector{Float64}(undef, nnz_estimate)
    k = 0 # Counter for current non-zero element

    b_phi_seconder = zeros(M^2)
    
    # Pre-calculate loop-invariant values
    prefactor = 1.0 - mutation_rate
    mut_term = -mutation_rate / numOfStrategies

    # Main loop: populate I, J, V vectors instead of modifying the matrix directly
    for i in 1:M
        Li = L_vector[i]
        for j in i+1:M
            Lj = L_vector[j]

            index1 = (i - 1) * M + j
            index2 = (j - 1) * M + i
            
            b_phi_seconder[index1] = mut_term
            # b_phi_seconder[index2] is already 0

            denom = Li + Lj
            
            # Add the first component for the value at (index1, index1)
            k += 1; I[k] = index1; J[k] = index1; V[k] = -1.0

            # Set values for (index2, index2) and (index2, index1)
            k += 1; I[k] = index2; J[k] = index2; V[k] = 1.0
            k += 1; I[k] = index2; J[k] = index1; V[k] = -1.0

            # Inner loop, continue adding elements to I, J, V
            inv_denom_prefactor = prefactor / denom
            for h in 1:M
                k += 1; I[k] = index1; J[k] = (h - 1) * M + j; V[k] = Lj * A[i, h] * inv_denom_prefactor
                k += 1; I[k] = index1; J[k] = (i - 1) * M + h; V[k] = Li * A[j, h] * inv_denom_prefactor
            end
        end
    end

    for i in 1:M
        index = (i - 1) * M + i
        b_phi_seconder[index] = 1.0
        k += 1; I[k] = index; J[k] = index; V[k] = 1.0
    end

    resize!.((I, J, V), k)

    coeff_matrix_phi_secondOrder = sparse(I, J, V, M^2, M^2)

    phi_secondOrder = idrs(coeff_matrix_phi_secondOrder, b_phi_seconder)
    phi_secondOrder_matrix = reshape(phi_secondOrder, (M, M))
    return phi_secondOrder_matrix
end


function multiple_calculate_phi_thirdOrder(A, num_nodes, offset, group_list, mutation_rate, numOfStrategies, phi_secondOrder_matrix)
    M = sum(group_list)
    M_sq = M^2
    M_cb = M^3

    # Efficiently build L_vector
    L_vector = Vector{Int}(undef, M)
    current_pos = 1
    for group_size in group_list
        range_end = current_pos + group_size - 1
        L_vector[current_pos:range_end] .= group_size
        current_pos = range_end + 1
    end

    # Initialize I, J, V vectors for COO format.
    nnz_estimate = round(Int, (M*(M-1)*(M-2)/6) * (12 + 3*M) + 3*M^2 + M)
    I = Vector{Int}(undef, nnz_estimate)
    J = Vector{Int}(undef, nnz_estimate)
    V = Vector{Float64}(undef, nnz_estimate)
    k_nnz = 0 # Counter for non-zero elements

    b_phi_thirdOrder = zeros(M_cb)

    # Pre-calculate invariants
    prefactor = 1.0 - mutation_rate
    mut_div_numstrat = mutation_rate / numOfStrategies


    #println("Filling I, J, V vectors:")
    # Main loop (i, j, k are all distinct)
    for i in 1:M
        Li = L_vector[i]
        inv_Li = 1.0 / Li
        for j in i+1:M
            Lj = L_vector[j]
            inv_Lj = 1.0 / Lj
            for k in j+1:M
                Lk = L_vector[k]
                inv_Lk = 1.0 / Lk

                idx_base = ((i-1)*M + j-1)*M + k

                b_phi_thirdOrder[idx_base] = -mut_div_numstrat * (phi_secondOrder_matrix[j, k] * inv_Li + phi_secondOrder_matrix[i, k] * inv_Lj + phi_secondOrder_matrix[i, j] * inv_Lk)

                k_nnz += 1; I[k_nnz] = idx_base; J[k_nnz] = idx_base; V[k_nnz] = -(inv_Li + inv_Lj + inv_Lk)
                
                # Symmetry equations
                permutations = [
                    ((i-1)*M + k-1)*M + j, # i,k,j
                    ((j-1)*M + i-1)*M + k, # j,i,k
                    ((j-1)*M + k-1)*M + i, # j,k,i
                    ((k-1)*M + i-1)*M + j, # k,i,j
                    ((k-1)*M + j-1)*M + i  # k,j,i
                ]
                for idx_p in permutations
                    k_nnz += 1; I[k_nnz] = idx_p; J[k_nnz] = idx_p; V[k_nnz] = 1.0
                    k_nnz += 1; I[k_nnz] = idx_p; J[k_nnz] = idx_base; V[k_nnz] = -1.0
                end

                # Inner h loop, leveraging sparsity of A
                prefactor_inv_Li = prefactor * inv_Li
                prefactor_inv_Lj = prefactor * inv_Lj
                prefactor_inv_Lk = prefactor * inv_Lk
                for h in 1:M
                    if A[i, h] != 0.0
                        k_nnz += 1; I[k_nnz] = idx_base; J[k_nnz] = ((h-1)*M + j-1)*M + k; V[k_nnz] = A[i, h] * prefactor_inv_Li
                    end
                    if A[j, h] != 0.0
                        k_nnz += 1; I[k_nnz] = idx_base; J[k_nnz] = ((i-1)*M + h-1)*M + k; V[k_nnz] = A[j, h] * prefactor_inv_Lj
                    end
                    if A[k, h] != 0.0
                        k_nnz += 1; I[k_nnz] = idx_base; J[k_nnz] = ((i-1)*M + j-1)*M + h; V[k_nnz] = A[k, h] * prefactor_inv_Lk
                    end
                end
            end
        end
    end

    # Loop for boundary conditions (two indices are the same)
    for i in 1:M
        for j in 1:M
            if i == j; continue; end
            val = phi_secondOrder_matrix[i, j]
            
            indices_to_set = [
                ((i-1)*M + i-1)*M + j, # (i,i,j)
                ((i-1)*M + j-1)*M + i, # (i,j,i)
                ((j-1)*M + i-1)*M + i, # (j,i,i) - This seems more symmetric than (i,j,j)
            ]
            
            for idx in unique(indices_to_set) # Use unique to handle cases like (i,i,j) == (j,i,i) if j<i
                 b_phi_thirdOrder[idx] = val
                 k_nnz += 1; I[k_nnz] = idx; J[k_nnz] = idx; V[k_nnz] = 1.0
            end
        end
    end

    # Boundary condition (all indices are the same)
    for i in 1:M
        index_diag = ((i-1)*M + i-1)*M + i
        b_phi_thirdOrder[index_diag] = 1.0
        k_nnz += 1; I[k_nnz] = index_diag; J[k_nnz] = index_diag; V[k_nnz] = 1.0
    end

    # Trim and build the sparse matrix
    resize!.((I, J, V), k_nnz)
    coeff_matrix_phi_thirdOrder = sparse(I, J, V, M_cb, M_cb)

    # Build initial guess vector x0 to accelerate convergence.
    # Initialize with random numbers in [0,1] as they are probabilities.
    x0 = rand(M_cb)
    for i in 1:M
        # Set exact values for known boundary conditions, overwriting random values.
        # Case i=j=k -> phi_iii = 1
        x0[((i-1)*M + i-1)*M + i] = 1.0
        for j in 1:M
            if i == j; continue; end
            # Case two indices are equal -> phi_iij = phi_ij
            val = phi_secondOrder_matrix[i, j]
            indices_to_set = [
                ((i-1)*M + i-1)*M + j, # (i,i,j)
                ((i-1)*M + j-1)*M + i, # (i,j,i)
                ((j-1)*M + i-1)*M + i, # (j,i,i)
            ]
            for idx in unique(indices_to_set)
                x0[idx] = val
            end
        end
    end
    #end

    # Solve and return (using x0 as the initial guess).
    # Use the in-place version idrs!, which takes the initial guess as the first argument 
    # and stores the result directly in x0.
    #println("Solving the linear system:")
    phi_thirdOrder = idrs!(x0, coeff_matrix_phi_thirdOrder, b_phi_thirdOrder)
    return phi_thirdOrder
end


function multiple_calculate_structure_coefficient(edge_list, num_nodes, group_list, group_data, numOfStrategies, mutation_rate, alpha=0)
    offset = calculate_offset(group_list)
    D = calculate_D(edge_list, num_nodes, group_list, group_data, offset)
    A = calculate_A(D)
    adj = edgeList2adjacencyMatrix(edge_list, num_nodes)
    W = adj2highorder_adj(adj, num_nodes, group_list, offset, group_data)
    M = sum(group_list)
    
    pi_mut = multiple_calculate_pi_mut(A, num_nodes, offset, group_list, mutation_rate)
    omega_matrix = multiple_calculate_omega(W, num_nodes, offset, group_list)
    
    #println("phi_secondOrder_matrix calculation:")
    phi_secondOrder_matrix = multiple_calculate_phi_secondOrder(A, num_nodes, offset, group_list, mutation_rate, numOfStrategies)
    
    #println("phi_thirdOrder calculation:")
    phi_thirdOrder = multiple_calculate_phi_thirdOrder(A, num_nodes, offset, group_list, mutation_rate, numOfStrategies, phi_secondOrder_matrix)
    
    # --- 2. Pre-computation for vectorization ---
    N = numOfStrategies
    mu = mutation_rate
    
    # Reshape third-order probabilities from a vector to a 3D tensor for easier indexing
    phi_3_tensor = reshape(phi_thirdOrder, (M, M, M))

    # Create helper vector to map a global individual index back to its node index
    node_of_member = Vector{Int}(undef, M)
    for i in 1:num_nodes
        node_of_member[(offset[i]+1):(offset[i]+group_list[i])] .= i
    end

    # Pre-calculate the sum of columns of A within each node
    A_node_sum = zeros(M, num_nodes)
    for k_node in 1:num_nodes
        node_members = (1:group_list[k_node]) .+ offset[k_node]
        A_node_sum[:, k_node] = sum(A[:, node_members], dims=2)
    end
    
    # Pre-calculate terms involving sums over omega_matrix
    omega_row_sum = sum(omega_matrix, dims=2)
    # S_o_p2_uw[u,w] = sum_x omega[w,x] * phi_2[u,x]
    S_o_p2_uw = phi_secondOrder_matrix * transpose(omega_matrix)
    # S_o_p2_ww[w] = sum_x omega[w,x] * phi_2[w,x]
    S_o_p2_ww = sum(omega_matrix .* phi_secondOrder_matrix, dims=2)

    # --- 3. Main calculation: replace O(M^8) with O(M^3) loops ---
    lambda1, lambda2, lambda3 = 0.0, 0.0, 0.0
    
    #println("lambda1, lambda2, lambda3 calculation (optimized):")
        # Define constants from the lambda expressions
        C_denom = mutation_rate * (N - 1) * (N - 2)
        mu_term = N * mu
        one_minus_mu = 1 - mu
        N_sq = N^2
        N_sq_1_mu = N_sq * one_minus_mu


        # Loop over u, v, w (global individual indices)
        @inbounds for w in 1:M
            w_node = node_of_member[w]
            
            # Pre-calculate sums over x (the innermost loop)
            # S_o_p3_uw[u] = sum_x omega[w,x] * phi_3[u,w,x]
            # S_o_p3_vw[v] = sum_x omega[w,x] * phi_3[v,w,x]
            S_o_p3_uw = phi_3_tensor[:, w, :] * omega_matrix[w, :]
            S_o_p3_vw = phi_3_tensor[:, w, :] * omega_matrix[w, :] # Same calculation, used for v index


            
            for v in 1:M
                v_node = node_of_member[v]
                
                # Phi_2 terms that depend on v and w
                phi2_vw = phi_secondOrder_matrix[v, w]
                S_o_p2_vw_val = S_o_p2_uw[v, w] # S_o_p2_uw is symmetric in usage for u and v
                
                @inbounds for u in 1:M 
                    A_uv = A[u, v]
                    if A_uv == 0.0 continue end

                    u_node = node_of_member[u]
                    
                    m_val = A_uv * ( (v_node == w_node) - A_node_sum[u, w_node] )
                    if m_val == 0.0 continue end

                    term_common = pi_mut[u] / (num_nodes * group_list[u_node]) * m_val
                    if term_common == 0.0 continue end
                    

                    phi2_uw = phi_secondOrder_matrix[u, w]
                    phi2_wx = S_o_p2_ww[w]
                    
                    l1_val = ( -N_sq * S_o_p3_uw[u] + N_sq_1_mu * S_o_p3_vw[v] 
                             + N * phi2_uw * omega_row_sum[w] + N * S_o_p2_uw[u, w] 
                             - N * one_minus_mu * phi2_vw * omega_row_sum[w] - N * one_minus_mu * S_o_p2_vw_val
                             + mu_term * phi2_wx - 2 * mu )
                    
                    l2_val = ( -N_sq * S_o_p3_uw[u] + N_sq_1_mu * S_o_p3_vw[v] 
                             + N * phi2_uw * omega_row_sum[w] + N * (N-1) * S_o_p2_uw[u, w] 
                             - N * one_minus_mu * phi2_vw * omega_row_sum[w] - N * (N-1) * one_minus_mu * S_o_p2_vw_val
                             + mu_term * phi2_wx - mu_term )

                    l3_val = ( 2*N_sq * S_o_p3_uw[u] - 2*N_sq_1_mu * S_o_p3_vw[v] 
                             - N_sq * phi2_uw * omega_row_sum[w] - N_sq * S_o_p2_uw[u, w] 
                             + N_sq_1_mu * phi2_vw * omega_row_sum[w] + N_sq_1_mu * S_o_p2_vw_val
                             - 2 * mu_term * phi2_wx + 2 * mu_term )

                    lambda1 += term_common * l1_val
                    lambda2 += term_common * l2_val
                    lambda3 += term_common * l3_val
                end
            end
        end

    scaling_factor = 1.0 / C_denom
    lambda1 *= scaling_factor
    lambda2 *= scaling_factor
    lambda3 *= scaling_factor

    return lambda1, lambda2, lambda3
end

function ultimatumGame_payoffMatrix(k)
    PM = zeros(k^2,k^2)
    for py in 0:k-1
        for qy in 0:k-1
            for pz in 0:k-1
                for qz in 0:k-1
                    if (py >= qz) && (pz >=qy)
                        PM[(py)*k+qy+1, (pz)*k+qz+1] = 1-py/(k-1) + pz/(k-1)
                    elseif (py >= qz) && (pz <qy)
                        PM[(py)*k+qy+1, (pz)*k+qz+1] = 1-py/(k-1)
                    elseif (py < qz) && (pz >=qy)
                        PM[(py)*k+qy+1, (pz)*k+qz+1] = pz/(k-1)
                    else
                        PM[(py)*k+qy+1, (pz)*k+qz+1] = 0
                    end
                end
            end
        end
    end
    return PM
end

function ultimatumGame_payoffProcessing(k)
    a_yStar_bar = zeros(k^2)
    a_Stary_bar = zeros(k^2)
    a_bar_1 = 0
    a_bar_2 = 0
    for py in 0:k-1
        for qy in 0:k-1
            a_yStar_bar[py*k+qy+1] = ( (k^2 - (k^2-k)*qy/(k-1))*(1+qy/(k-1))/2 + (k + (k^2-k)*py/(k-1))*(1-py/(k-1)))/k^2
            a_Stary_bar[py*k+qy+1] = ( (k^2 - (k^2-k)*qy/(k-1))*(1-qy/(k-1))/2 + (k + (k^2-k)*py/(k-1))*(py/(k-1)))/k^2
        end
        a_bar_1 += (k^2/2+(3/2*k-k^2)*py/(k-1)+(3/2*k^2-3/2*k)*(py/(k-1))^2)/k^3
        a_bar_2 += (k^2/2 + k -(3/2*k-k^2)*py/(k-1)-(3/2*k^2-3/2*k)*(py/(k-1))^2)/k^3
    end
    return a_yStar_bar, a_Stary_bar,a_bar_1 ,a_bar_2
end

function ultimatumGame_pq_frequency_discrete(lambda1, lambda2, lambda3, k, selection_intensity)
    p = 1/2 + selection_intensity * ( lambda1*(1/12 + 1/(12*k)) + lambda2*(-5/12 +7/(12*k) + (k-3)*(2*k-1)/(6*k*(k-1))) + lambda3*(1//(12*k) - 1/(6*(k-1))) )
    q = 1/2 + selection_intensity * ( lambda1*(-1/12-1/(12*k))  + lambda2*(1/(6*(k-1)) - 1/(12*k))                      + lambda3*(1/(12*(k-1)) - 1/(12*k) - 1/24) )
    return p,q
end

function ultimatumGame_pq_frequency_continual(lambda1, lambda2, lambda3, selection_intensity)
    p = 1/2 + selection_intensity * ( lambda1*1/12     - lambda2*1/12)
    q = 1/2 + selection_intensity * ( lambda1*(-1/12)  - lambda3*(-1/24) )
    return p,q
end

function multiple_donationGame_K1K2(edge_list, num_nodes, group_list, group_data, mutation_rate, alpha)
    offset = calculate_offset(group_list)
    D = calculate_D(edge_list, num_nodes, group_list, group_data, offset)
    A = calculate_A(D)
    adj = edgeList2adjacencyMatrix(edge_list, num_nodes)
    W = adj2highorder_adj(adj, num_nodes, group_list, offset, group_data)
    M = sum(group_list)
    omega_matrix = multiple_calculate_omega(W, num_nodes, offset, group_list)

    omega_sum = sum(omega_matrix, dims=2)
    omega_sum_columnMatrix = repeat(reshape(omega_sum, (1,M)), inner=(M, 1))

    L_vector = vcat([repeat([x], x) for x in group_list]...)
    
    pi_mut = multiple_calculate_pi_mut(A, num_nodes, offset, group_list, mutation_rate)

    pi_over_L = reshape(pi_mut./L_vector, (M,1))
    pi_over_L_rowMatrix = repeat(pi_over_L, inner=(1, M))

    A_imjn_times_pi_im_over_L_im = pi_over_L_rowMatrix .* A

    delta_jn_kl_matrix = zeros(M,M)
    for i in 1:num_nodes
        m_offset = offset[i]
        group_size = group_list[i]
        delta_jn_kl_matrix[m_offset+1:m_offset+group_size, m_offset+1:m_offset+group_size] .= 1
    end

    A_sum_l = zeros(M, M)
    for i in 1:M
        for k in 1:num_nodes
            sum_A_im_kl_for_l = sum(A[i, offset[k]+1: (offset[k]+group_list[k])])
            A_sum_l[i, offset[k]+1: (offset[k]+group_list[k])] .= sum_A_im_kl_for_l
        end
    end

    K11 = 0
    K12 = 0
    K13 = 0
    K14 = 0
    K15 = 0
    K16 = 0

    K21 = 0
    K22 = 0
    K23 = 0
    K24 = 0
    K25 = 0
    K26 = 0

    phi_secondOrder_matrix = multiple_calculate_phi_secondOrder(A, num_nodes, offset, group_list, mutation_rate, 2)

    K11_temp = omega_sum_columnMatrix .* phi_secondOrder_matrix .* delta_jn_kl_matrix
    K11 = sum(A_imjn_times_pi_im_over_L_im * K11_temp) *2 * (1-mutation_rate) /(mutation_rate * num_nodes)

    K12_temp1 = omega_sum_columnMatrix .* phi_secondOrder_matrix
    K12_temp2 = A_imjn_times_pi_im_over_L_im * delta_jn_kl_matrix
    K12 = sum(K12_temp1 .* K12_temp2) * 2 / (num_nodes * mutation_rate)

    K13_temp = omega_sum_columnMatrix .* delta_jn_kl_matrix
    K13 = sum(A_imjn_times_pi_im_over_L_im * K13_temp) / num_nodes

    K14_temp = A_imjn_times_pi_im_over_L_im * phi_secondOrder_matrix
    K14 = sum(K14_temp .* omega_sum_columnMatrix .* A_sum_l) *2 *(1-mutation_rate) / (num_nodes*mutation_rate)

    K15_temp = A_sum_l .* phi_secondOrder_matrix .* omega_sum_columnMatrix
    K15 = sum(sum(K15_temp, dims=2) .* sum(A_imjn_times_pi_im_over_L_im, dims=2)) *2/(num_nodes *mutation_rate)

    K16_temp = A_sum_l .*omega_sum_columnMatrix
    K16 = sum( sum(K16_temp, dims=2) .* sum(A_imjn_times_pi_im_over_L_im, dims=2)) /num_nodes

    K21_temp1 = delta_jn_kl_matrix * omega_matrix
    K21_temp2 = K21_temp1 .* phi_secondOrder_matrix
    K21 = sum(A_imjn_times_pi_im_over_L_im * K21_temp2) * 2 * (1-mutation_rate)/ (num_nodes*mutation_rate)

    K22_temp1 = delta_jn_kl_matrix * omega_matrix
    K22_temp2 = A_imjn_times_pi_im_over_L_im * K22_temp1
    K22 = sum(K22_temp2 .* phi_secondOrder_matrix) *2 / (num_nodes * mutation_rate)

    K23 = sum(A_imjn_times_pi_im_over_L_im * delta_jn_kl_matrix * omega_matrix)/num_nodes

    K24_temp1 = omega_matrix * phi_secondOrder_matrix
    K24_temp2 = A_sum_l * K24_temp1
    K24 = sum(K24_temp2 .* A_imjn_times_pi_im_over_L_im) *2 * (1-mutation_rate) / (num_nodes * mutation_rate)


    K25_temp1 = omega_matrix * phi_secondOrder_matrix
    K25_temp2 = A_sum_l .* transpose(K25_temp1)
    K25 = sum(sum(A_imjn_times_pi_im_over_L_im, dims=2) .* sum(K25_temp2, dims=2)) *2/(num_nodes*mutation_rate)

    K26_temp1 = A_sum_l * omega_matrix
    K26 = sum(sum(A_imjn_times_pi_im_over_L_im, dims=2) .* sum(K26_temp1, dims=2)) /num_nodes


    K1 = K11 - K12 + K13 - K14 + K15 - K16
    K2 = K21 - K22 + K23 - K24 + K25 - K26

    return K1, K2
end 


#MS-BD
function multiple_calculate_structure_coefficient_BD(edge_list, num_nodes, group_list, group_data, numOfStrategies, mutation_rate, alpha=0)
    offset = calculate_offset(group_list)
    D = calculate_D(edge_list, num_nodes, group_list, group_data, offset)
    adj = edgeList2adjacencyMatrix(edge_list, num_nodes)
    W = adj2highorder_adj(adj, num_nodes, group_list, offset, group_data)
    M = sum(group_list)
    E = calculate_E(D)
    I = Diagonal(ones(num_nodes)) 

    omega_matrix = multiple_calculate_omega(W, num_nodes, offset, group_list)
    Ndeath_rate_vector = multiple_calculate_Ndeath_rate_vector(M, num_nodes, offset, group_list, E)
    pi_mut = multiple_calculate_pi_mut_BD(E, num_nodes, offset, group_list, mutation_rate, Ndeath_rate_vector)
    phi_secondOrder_matrix = multiple_calculate_phi_secondOrder_BD(E, num_nodes, offset, group_list, mutation_rate, numOfStrategies, Ndeath_rate_vector)
    phi_thirdOrder = multiple_calculate_phi_thirdOrder_BD(E, num_nodes, offset, group_list, mutation_rate, numOfStrategies, phi_secondOrder_matrix, Ndeath_rate_vector)


    lambda1 = 0
    lambda2 = 0
    lambda3 = 0

    for i in 1:num_nodes
        for m in 1:group_list[i]
            for j in 1:num_nodes
                Lj = group_list[j]
                for n in 1:group_list[j]
                    for k in 1:num_nodes
                        for l in 1:group_list[k]
                            for h in 1:num_nodes
                                for s in 1:group_list[h]
                                    lambda1 += 1/mutation_rate * pi_mut[m+offset[i]] * ( 1/(num_nodes*Lj) * E[n + offset[j], m+offset[i]] * (I[j,k] - 1/num_nodes)) * omega_matrix[l + offset[k], s+offset[h]] * ( -numOfStrategies^2 * phi_thirdOrder[((m + offset[i] -1) * M + (l + offset[k] ) -1)*M + (s + offset[h])] + numOfStrategies^2 * (1-mutation_rate)* phi_thirdOrder[((n + offset[j] -1) * M + (l + offset[k] ) -1)*M + (s + offset[h])] + numOfStrategies*phi_secondOrder_matrix[m+offset[i], l+offset[k]] +  numOfStrategies*phi_secondOrder_matrix[m+offset[i], s+offset[h]] - numOfStrategies* (1-mutation_rate)*phi_secondOrder_matrix[n+offset[j], l+offset[k]] - numOfStrategies* (1-mutation_rate)*phi_secondOrder_matrix[n+offset[j], s+offset[h]] + numOfStrategies*mutation_rate*phi_secondOrder_matrix[l+offset[k], s+offset[h]] -2*mutation_rate )/( (numOfStrategies-1)*(numOfStrategies-2))
                                    lambda2 += 1/mutation_rate * pi_mut[m+offset[i]] * ( 1/(num_nodes*Lj) * E[n + offset[j], m+offset[i]] * (I[j,k] - 1/num_nodes)) * omega_matrix[l + offset[k], s+offset[h]] * ( -numOfStrategies^2 * phi_thirdOrder[((m + offset[i] -1) * M + (l + offset[k] ) -1)*M + (s + offset[h])] + numOfStrategies^2 * (1-mutation_rate)* phi_thirdOrder[((n + offset[j] -1) * M + (l + offset[k] ) -1)*M + (s + offset[h])] + numOfStrategies*phi_secondOrder_matrix[m+offset[i], l+offset[k]] +  numOfStrategies*(numOfStrategies-1)*phi_secondOrder_matrix[m+offset[i], s+offset[h]] - numOfStrategies* (1-mutation_rate)*phi_secondOrder_matrix[n+offset[j], l+offset[k]] - numOfStrategies*(numOfStrategies-1)* (1-mutation_rate)*phi_secondOrder_matrix[n+offset[j], s+offset[h]] + numOfStrategies*mutation_rate*phi_secondOrder_matrix[l+offset[k], s+offset[h]] -numOfStrategies*mutation_rate )/( (numOfStrategies-1)*(numOfStrategies-2))
                                    lambda3 += 1/mutation_rate * pi_mut[m+offset[i]] * ( 1/(num_nodes*Lj) * E[n + offset[j], m+offset[i]] * (I[j,k] - 1/num_nodes)) * omega_matrix[l + offset[k], s+offset[h]] * ( 2*numOfStrategies^2 * phi_thirdOrder[((m + offset[i] -1) * M + (l + offset[k] ) -1)*M + (s + offset[h])] - 2*numOfStrategies^2 * (1-mutation_rate)* phi_thirdOrder[((n + offset[j] -1) * M + (l + offset[k] ) -1)*M + (s + offset[h])] - numOfStrategies^2*phi_secondOrder_matrix[m+offset[i], l+offset[k]] -  numOfStrategies^2*phi_secondOrder_matrix[m+offset[i], s+offset[h]] + numOfStrategies^2 * (1-mutation_rate)*phi_secondOrder_matrix[n+offset[j], l+offset[k]] + numOfStrategies^2* (1-mutation_rate)*phi_secondOrder_matrix[n+offset[j], s+offset[h]] - 2*numOfStrategies*mutation_rate*phi_secondOrder_matrix[l+offset[k], s+offset[h]] +2*numOfStrategies*mutation_rate  )/( (numOfStrategies-1)*(numOfStrategies-2))
                                end
                            end
                        end
                    end
                end
            end
        end
    end

    return lambda1, lambda2, lambda3
end

function multiple_donationGame_K1K2_BD(edge_list, num_nodes, group_list, group_data, mutation_rate, alpha=0)
    offset = calculate_offset(group_list)
    D = calculate_D(edge_list, num_nodes, group_list, group_data, offset)
    adj = edgeList2adjacencyMatrix(edge_list, num_nodes)
    W = adj2highorder_adj(adj, num_nodes, group_list, offset, group_data)
    M = sum(group_list)
    E = calculate_E(D)
    I = Diagonal(ones(num_nodes)) 
    omega_matrix = multiple_calculate_omega(W, num_nodes, offset, group_list)

    K1 = 0
    K2 = 0

    if mutation_rate > 0
        Ndeath_rate_vector = multiple_calculate_Ndeath_rate_vector(M, num_nodes, offset, group_list, E)
        pi_mut = multiple_calculate_pi_mut_BD(E, num_nodes, offset, group_list, mutation_rate, Ndeath_rate_vector)
        phi_secondOrder_matrix = multiple_calculate_phi_secondOrder_BD(E, num_nodes, offset, group_list, mutation_rate, 2, Ndeath_rate_vector)
        #println("phi second order:$(phi_secondOrder_matrix)")
        for i in 1:num_nodes
            for m in 1:group_list[i]
                for j in 1:num_nodes
                    Lj = group_list[j]
                    for n in 1:group_list[j]
                        for k in 1:num_nodes
                            for l in 1:group_list[k]
                                for h in 1:num_nodes
                                    for s in 1:group_list[h]
                                        K1 += 1/mutation_rate * pi_mut[m + offset[i]] * ( 1/(num_nodes*Lj) * E[n + offset[j], m+offset[i]] * (I[j,k] - 1/num_nodes)) * omega_matrix[l + offset[k], s + offset[h]] * ( 2*(1-mutation_rate)*phi_secondOrder_matrix[n+offset[j], l + offset[k]] - 2 *phi_secondOrder_matrix[m + offset[i], l + offset[k]] + mutation_rate )
                                        K2 += 1/mutation_rate * pi_mut[m + offset[i]] * ( 1/(num_nodes*Lj) * E[n + offset[j], m+offset[i]] * (I[j,k] - 1/num_nodes)) * omega_matrix[l + offset[k], s + offset[h]] * ( 2*(1-mutation_rate)*phi_secondOrder_matrix[n+offset[j], s + offset[h]] - 2 *phi_secondOrder_matrix[m + offset[i], s + offset[h]] + mutation_rate )
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end

    else
        println("Wrong mutation rate")
        return 0
    end

    return K1, K2

end


function multiple_donationGame_K1K2_BD_new(edge_list, num_nodes, group_list, group_data, mutation_rate, alpha=0)
    offset = calculate_offset(group_list)
    D = calculate_D(edge_list, num_nodes, group_list, group_data, offset)
    adj = edgeList2adjacencyMatrix(edge_list, num_nodes)
    W = adj2highorder_adj(adj, num_nodes, group_list, offset, group_data)
    M = sum(group_list)
    E = calculate_E(D)
    I = Diagonal(ones(num_nodes)) 
    omega_matrix = multiple_calculate_omega(W, num_nodes, offset, group_list)

    K1 = 0
    K2 = 0

    Ndeath_rate_vector = multiple_calculate_Ndeath_rate_vector(M, num_nodes, offset, group_list, E)
    pi_mut = multiple_calculate_pi_mut_BD(E, num_nodes, offset, group_list, mutation_rate, Ndeath_rate_vector)
    phi_secondOrder_matrix = multiple_calculate_phi_secondOrder_BD(E, num_nodes, offset, group_list, mutation_rate, 2, Ndeath_rate_vector)

    delta_matrix = zeros(M,M)
    for i in 1:num_nodes
        m_offset = offset[i]
        group_size = group_list[i]
        delta_matrix[m_offset+1:m_offset+group_size, m_offset+1:m_offset+group_size] .= 1
    end

    L_vector = vcat([repeat([x], x) for x in group_list]...)
    pi_matrix = repeat(pi_mut', length(pi_mut), 1) 
    L_matrix = repeat(L_vector, 1, M)  
    M1 = pi_matrix ./L_matrix .* E

    K11 = sum(transpose(M1) * (delta_matrix.*phi_secondOrder_matrix) * omega_matrix ) *2 *(1-mutation_rate) / (num_nodes*mutation_rate)
    K12 = sum(transpose(M1) * (phi_secondOrder_matrix * omega_matrix)) *2 *(1-mutation_rate) / (num_nodes^2*mutation_rate)
    K13 = sum(((transpose(M1) * delta_matrix) .* phi_secondOrder_matrix) * omega_matrix ) *2/ (num_nodes*mutation_rate)
    K14 = sum(M1 * phi_secondOrder_matrix * omega_matrix) *2 / (num_nodes^2*mutation_rate)
    K15 = sum(transpose(M1) * delta_matrix * omega_matrix) /num_nodes
    K16 = sum(M1)*sum(omega_matrix)/num_nodes^2

    K21 = sum( transpose(M1) * ((delta_matrix * omega_matrix) .* phi_secondOrder_matrix)) *2 *(1-mutation_rate) / (num_nodes*mutation_rate)
    K22 = sum(omega_matrix * transpose(phi_secondOrder_matrix) * M1) *2 *(1-mutation_rate) / (num_nodes^2*mutation_rate)
    K23 = sum( (M1 * phi_secondOrder_matrix) .* (delta_matrix * omega_matrix) )*2/ (num_nodes*mutation_rate)
    K24 = sum(M1 * phi_secondOrder_matrix * transpose(omega_matrix)) *2 / (num_nodes^2*mutation_rate)
    K25 = sum(transpose(M1)  * delta_matrix * omega_matrix) /num_nodes
    K26 = K16

    K1 = K11 - K12 - K13 + K14 + K15 - K16
    K2 = K21 - K22 - K23 + K24 + K25 - K26

    return K1, K2

end


function multiple_calculate_pi_mut_BD(E, num_nodes, offset, group_list, mutation_rate, Ndeath_rate_vector)
    M = sum(group_list)
    coeff_matrix = zeros(M,M)
    b_vector = zeros(M,1)

    for i in 1:num_nodes
        Li = group_list[i]
        for m in 1:group_list[i]
            b_vector[m+offset[i]] = - mutation_rate * num_nodes / M 
            coeff_matrix[m+offset[i], m+offset[i]] = - Ndeath_rate_vector[m+offset[i]]
            for j in 1:num_nodes
                for n in 1:group_list[j]
                    coeff_matrix[m+offset[i], n+offset[j]] += (1-mutation_rate)/Li * E[m+offset[i], n+offset[j]] 
                end
            end
        end
    end
   
    #pi_mut_vector = coeff_matrix \ b_vector
    pi_mut_vector = idrs(coeff_matrix,b_vector) 
    return pi_mut_vector
end

function multiple_calculate_Ndeath_rate_vector(M, num_nodes, offset, group_list, E)
    L_vector = vcat([repeat([x], x) for x in group_list]...)
    num_columns = length(M)  
    L_matrix = repeat(L_vector, 1, num_columns)
    E_over_lk = E ./ L_matrix
    Ndeath_rate_vector = sum(E_over_lk, dims=1)  

    return Ndeath_rate_vector
end

function multiple_calculate_Ndeath_rate_vector_ori(M, num_nodes, offset, group_list, E)
    Ndeath_rate_vector = zeros(M)

    for i in 1:num_nodes
        for m in 1:group_list[i]
            for k in 1:num_nodes
                for l in 1:group_list[k]
                    Ndeath_rate_vector[m+offset[i]] += E[l+offset[k], m+offset[i]] / group_list[k]
                end
            end
        end
    end

    return Ndeath_rate_vector
end

function multiple_calculate_phi_secondOrder_BD(E, num_nodes, offset, group_list, mutation_rate, numOfStrategies, Ndeath_rate_vector)
    M = sum(group_list)
    coeff_matrix_phi_secondOrder = spzeros(M^2, M^2)
    b_phi_seconder = zeros(M^2)

    for i in 1:num_nodes
        for m in 1:group_list[i]
            for j in 1:num_nodes
                for n in 1:group_list[j]
                    index = (m + offset[i] -1) * M + (n + offset[j] )

                    if (i==j)&&(m==n)
                        coeff_matrix_phi_secondOrder[index,index] = 1
                        b_phi_seconder[index] = 1
                    else
                        coeff_matrix_phi_secondOrder[index,index] = -1
                        b_phi_seconder[index] = -mutation_rate / numOfStrategies

                        for h in 1:num_nodes
                            Lh = group_list[h]
                            for s in 1:group_list[h]
                                coeff_matrix_phi_secondOrder[index, (s + offset[h] -1) * M + (n + offset[j] )] += E[s+offset[h], m+offset[i]]/Lh * (1-mutation_rate)/ ( Ndeath_rate_vector[m+offset[i]] + Ndeath_rate_vector[n+offset[j]] )
                                coeff_matrix_phi_secondOrder[index, (m + offset[i] -1) * M + (s + offset[h] )] += E[s+offset[h], n+offset[j]]/Lh * (1-mutation_rate)/ ( Ndeath_rate_vector[m+offset[i]] + Ndeath_rate_vector[n+offset[j]] )
                            end
                        end
                    end
                end
            end
        end
    end

    #phi_secondOrder = coeff_matrix_phi_secondOrder \ b_phi_seconder
    phi_secondOrder = idrs(coeff_matrix_phi_secondOrder, b_phi_seconder)
    phi_secondOrder_matrix = reshape(phi_secondOrder, (M, M))
    return phi_secondOrder_matrix
end

function multiple_calculate_phi_thirdOrder_BD(E, num_nodes, offset, group_list, mutation_rate, numOfStrategies, phi_secondOrder_matrix, Ndeath_rate_vector)
    M = sum(group_list)
    coeff_matrix_phi_thirdOrder = spzeros(M^3, M^3)
    b_phi_thirdOrder = zeros(M^3)
    L_vector = vcat([repeat([x], x) for x in group_list]...)
    # 生成 Ndeath_rate_vector_resize
    #Ndeath_rate_vector_resize = vcat([repeat([x], n) for (x, n) in zip(Ndeath_rate_vector, group_list)]...)
    Ndeath_rate_vector_resize = Ndeath_rate_vector
    for i in 1:M
        for j in i+1:M
            for k in j+1:M
                # im - min jn - medium kl - max
                index1 = ((i -1) * M + (j ) -1)*M + (k) # i<j<k
                index2 = ((i -1) * M + (k ) -1)*M + (j) # i<k<j
                index3 = ((j -1) * M + (i ) -1)*M + (k) # j<i<K
                index4 = ((k -1) * M + (i ) -1)*M + (j) # j<k<i
                index5 = ((j -1) * M + (k ) -1)*M + (i) # k<i<j
                index6 = ((k -1) * M + (j ) -1)*M + (i) # k<j<i
                coeff_matrix_phi_thirdOrder[index1, index1] = -1
                coeff_matrix_phi_thirdOrder[index2, index2] = 1
                coeff_matrix_phi_thirdOrder[index3, index3] = 1
                coeff_matrix_phi_thirdOrder[index4, index4] = 1
                coeff_matrix_phi_thirdOrder[index5, index5] = 1
                coeff_matrix_phi_thirdOrder[index6, index6] = 1
                #b_phi_thirdOrder[index1] = -mutation_rate * ( phi_secondOrder_matrix[ j, k ]/(numOfStrategies * Li) + phi_secondOrder_matrix[ i, k ]/(numOfStrategies * Lj) + phi_secondOrder_matrix[i, j ]/(numOfStrategies * Lk))
                b_phi_thirdOrder[index1] = -mutation_rate/(numOfStrategies * (Ndeath_rate_vector_resize[i] + Ndeath_rate_vector_resize[j] + Ndeath_rate_vector_resize[k])) * ( Ndeath_rate_vector_resize[i]*phi_secondOrder_matrix[j, k] + Ndeath_rate_vector_resize[j]*phi_secondOrder_matrix[i, k] + Ndeath_rate_vector_resize[k]*phi_secondOrder_matrix[i, j] )
                b_phi_thirdOrder[index2] = 0
                b_phi_thirdOrder[index3] = 0
                b_phi_thirdOrder[index4] = 0
                b_phi_thirdOrder[index5] = 0
                b_phi_thirdOrder[index6] = 0
                coeff_matrix_phi_thirdOrder[index2, index1] = -1 
                coeff_matrix_phi_thirdOrder[index3, index1] = -1 
                coeff_matrix_phi_thirdOrder[index4, index1] = -1 
                coeff_matrix_phi_thirdOrder[index5, index1] = -1 
                coeff_matrix_phi_thirdOrder[index6, index1] = -1 

                for h in 1:M
                    Lh = L_vector[h]
                    coeff_matrix_phi_thirdOrder[index1, ((h -1) * M + (j ) -1)*M + (k)] += (1-mutation_rate) / (Lh * (Ndeath_rate_vector_resize[i] + Ndeath_rate_vector_resize[j] + Ndeath_rate_vector_resize[k])) * E[h,i]
                    coeff_matrix_phi_thirdOrder[index1, ((i -1) * M + (h ) -1)*M + (k)] += (1-mutation_rate) / (Lh * (Ndeath_rate_vector_resize[i] + Ndeath_rate_vector_resize[j] + Ndeath_rate_vector_resize[k])) * E[h,j]
                    coeff_matrix_phi_thirdOrder[index1, ((i -1) * M + (j ) -1)*M + (h)] += (1-mutation_rate) / (Lh * (Ndeath_rate_vector_resize[i] + Ndeath_rate_vector_resize[j] + Ndeath_rate_vector_resize[k])) * E[h,k]
                end
            end
        end
    end

    for i in 1:M
        for j in 1:M
            index1 = ((i -1) * M + (i ) -1)*M + (j) # i=j!=k
            index2 = ((i -1) * M + (j ) -1)*M + (i) # i=k!=j
            index3 = ((i -1) * M + (j ) -1)*M + (j) # k=j!=i
            coeff_matrix_phi_thirdOrder[index1, index1] = 1
            coeff_matrix_phi_thirdOrder[index2, index2] = 1
            coeff_matrix_phi_thirdOrder[index3, index3] = 1
            b_phi_thirdOrder[index1] = phi_secondOrder_matrix[ i, j ]
            b_phi_thirdOrder[index2] = phi_secondOrder_matrix[ i, j ]
            b_phi_thirdOrder[index3] = phi_secondOrder_matrix[ i, j ]
        end
        index_diag = ((i -1) * M + (i ) -1)*M + (i)
        b_phi_thirdOrder[index_diag] = 1
    end
    

    #phi_thirdOrder = coeff_matrix_phi_thirdOrder \ b_phi_thirdOrder
    phi_thirdOrder = idrs(coeff_matrix_phi_thirdOrder , b_phi_thirdOrder)
    return phi_thirdOrder
end


#MS-PC
function multiple_calculate_structure_coefficient_PC(edge_list, num_nodes, group_list, group_data, numOfStrategies, mutation_rate, alpha=0)
    offset = calculate_offset(group_list)
    D = calculate_D(edge_list, num_nodes, group_list, group_data, offset)
    A = calculate_A(D)
    adj = edgeList2adjacencyMatrix(edge_list, num_nodes)
    W = adj2highorder_adj(adj, num_nodes, group_list, offset, group_data)
    M = sum(group_list)

    pi_mut = multiple_calculate_pi_mut_PC(A, num_nodes, offset, group_list, mutation_rate)
    omega_matrix = multiple_calculate_omega(W, num_nodes, offset, group_list)
    phi_secondOrder_matrix = multiple_calculate_phi_secondOrder_PC(A, num_nodes, offset, group_list, mutation_rate, numOfStrategies)
    phi_thirdOrder = multiple_calculate_phi_thirdOrder_PC(A, num_nodes, offset, group_list, mutation_rate, numOfStrategies, phi_secondOrder_matrix)
    I = Diagonal(ones(num_nodes)) 

    #println("pi mut:$(pi_mut)\nphi second:$(phi_secondOrder_matrix)")

    lambda1 = 0
    lambda2 = 0
    lambda3 = 0

    for i in 1:num_nodes
        for m in 1:group_list[i]
            for j in 1:num_nodes
                for n in 1:group_list[j]
                    for k in 1:num_nodes
                        for l in 1:group_list[k]
                            for h in 1:num_nodes
                                for s in 1:group_list[h]
                                    lambda1 += 1/mutation_rate * pi_mut[m+offset[i]] * multiple_calculate_m_coefficient_PC(A, I, num_nodes, offset, group_list, alpha, k,l,j,n,i,m) * omega_matrix[l + offset[k], s+offset[h]] * ( -numOfStrategies^2 * phi_thirdOrder[((m + offset[i] -1) * M + (l + offset[k] ) -1)*M + (s + offset[h])] + numOfStrategies^2 * (1-mutation_rate)* phi_thirdOrder[((n + offset[j] -1) * M + (l + offset[k] ) -1)*M + (s + offset[h])] + numOfStrategies*phi_secondOrder_matrix[m+offset[i], l+offset[k]] +  numOfStrategies*phi_secondOrder_matrix[m+offset[i], s+offset[h]] - numOfStrategies* (1-mutation_rate)*phi_secondOrder_matrix[n+offset[j], l+offset[k]] - numOfStrategies* (1-mutation_rate)*phi_secondOrder_matrix[n+offset[j], s+offset[h]] + numOfStrategies*mutation_rate*phi_secondOrder_matrix[l+offset[k], s+offset[h]] -2*mutation_rate )/( (numOfStrategies-1)*(numOfStrategies-2))
                                    lambda2 += 1/mutation_rate * pi_mut[m+offset[i]] * multiple_calculate_m_coefficient_PC(A, I, num_nodes, offset, group_list, alpha, k,l,j,n,i,m) * omega_matrix[l + offset[k], s+offset[h]] * ( -numOfStrategies^2 * phi_thirdOrder[((m + offset[i] -1) * M + (l + offset[k] ) -1)*M + (s + offset[h])] + numOfStrategies^2 * (1-mutation_rate)* phi_thirdOrder[((n + offset[j] -1) * M + (l + offset[k] ) -1)*M + (s + offset[h])] + numOfStrategies*phi_secondOrder_matrix[m+offset[i], l+offset[k]] +  numOfStrategies*(numOfStrategies-1)*phi_secondOrder_matrix[m+offset[i], s+offset[h]] - numOfStrategies* (1-mutation_rate)*phi_secondOrder_matrix[n+offset[j], l+offset[k]] - numOfStrategies*(numOfStrategies-1)* (1-mutation_rate)*phi_secondOrder_matrix[n+offset[j], s+offset[h]] + numOfStrategies*mutation_rate*phi_secondOrder_matrix[l+offset[k], s+offset[h]] -numOfStrategies*mutation_rate )/( (numOfStrategies-1)*(numOfStrategies-2))
                                    lambda3 += 1/mutation_rate * pi_mut[m+offset[i]] * multiple_calculate_m_coefficient_PC(A, I, num_nodes, offset, group_list, alpha, k,l,j,n,i,m) * omega_matrix[l + offset[k], s+offset[h]] * ( 2*numOfStrategies^2 * phi_thirdOrder[((m + offset[i] -1) * M + (l + offset[k] ) -1)*M + (s + offset[h])] - 2*numOfStrategies^2 * (1-mutation_rate)* phi_thirdOrder[((n + offset[j] -1) * M + (l + offset[k] ) -1)*M + (s + offset[h])] - numOfStrategies^2*phi_secondOrder_matrix[m+offset[i], l+offset[k]] -  numOfStrategies^2*phi_secondOrder_matrix[m+offset[i], s+offset[h]] + numOfStrategies^2 * (1-mutation_rate)*phi_secondOrder_matrix[n+offset[j], l+offset[k]] + numOfStrategies^2* (1-mutation_rate)*phi_secondOrder_matrix[n+offset[j], s+offset[h]] - 2*numOfStrategies*mutation_rate*phi_secondOrder_matrix[l+offset[k], s+offset[h]] +2*numOfStrategies*mutation_rate  )/( (numOfStrategies-1)*(numOfStrategies-2))
                                end
                            end
                        end
                    end
                end
            end
        end
    end

    return lambda1, lambda2, lambda3
end

function multiple_calculate_pi_mut_PC(A, num_nodes, offset, group_list, mutation_rate)
    M = sum(group_list)
    coeff_matrix = zeros(M,M)
    b_vector = zeros(M,1)

    for i in 1:num_nodes
        for m in 1:group_list[i]
            b_vector[m+offset[i]] = -2 * mutation_rate * group_list[i] * num_nodes / (M * (1+mutation_rate)) 
            coeff_matrix[m+offset[i], m+offset[i]] = -1
            for j in 1:num_nodes
                for n in 1:group_list[j]
                    coeff_matrix[m+offset[i], n+offset[j]] += (1-mutation_rate)/(1+mutation_rate) * A[n+offset[j], m+offset[i]] * group_list[i] / group_list[j]
                end
            end
        end
    end
   
    #pi_mut_vector = coeff_matrix \ b_vector
    pi_mut_vector = idrs(coeff_matrix,b_vector) 
    
    return pi_mut_vector
end



function multiple_calculate_phi_secondOrder_PC(A, num_nodes, offset, group_list, mutation_rate, numOfStrategies)
    M = sum(group_list)
    #coeff_matrix = zeros(M^2, M^2)
    coeff_matrix_phi_secondOrder = spzeros(M^2, M^2)
    b_phi_seconder = zeros(M^2)

    for i in 1:num_nodes
        Li = group_list[i]
        for m in 1:group_list[i]
            for j in 1:num_nodes
                Lj = group_list[j]
                for n in 1:group_list[j]
                    index = (m + offset[i] -1) * M + (n + offset[j] )

                    if (i==j)&&(m==n)
                        coeff_matrix_phi_secondOrder[index,index] = 1
                        b_phi_seconder[index] = 1
                    else
                        coeff_matrix_phi_secondOrder[index,index] = -1
                        b_phi_seconder[index] = -2*mutation_rate / (numOfStrategies * (1+mutation_rate))

                        for h in 1:num_nodes
                            for s in 1:group_list[h]
                                coeff_matrix_phi_secondOrder[index, (s + offset[h] -1) * M + (n + offset[j] )] += Lj * A[m+offset[i], s+offset[h]] * (1-mutation_rate)/((Li + Lj) * (1+mutation_rate))
                                coeff_matrix_phi_secondOrder[index, (m + offset[i] -1) * M + (s + offset[h] )] += Li * A[n+offset[j], s+offset[h]] * (1-mutation_rate)/((Li + Lj) * (1+mutation_rate))
                            end
                        end
                    end
                end
            end
        end
    end

    #phi_secondOrder = coeff_matrix_phi_secondOrder \ b_phi_seconder
    phi_secondOrder = idrs(coeff_matrix_phi_secondOrder, b_phi_seconder)

    phi_secondOrder_matrix = reshape(phi_secondOrder, (M, M))
    return phi_secondOrder_matrix
end

function multiple_calculate_phi_thirdOrder_PC(A, num_nodes, offset, group_list, mutation_rate, numOfStrategies, phi_secondOrder_matrix)
    M = sum(group_list)
    coeff_matrix_phi_thirdOrder = spzeros(M^3, M^3)
    b_phi_thirdOrder = zeros(M^3)
    L_vector = vcat([repeat([x], x) for x in group_list]...)

    for i in 1:M
        Li = L_vector[i]
        for j in i+1:M
            Lj = L_vector[j]
            for k in j+1:M
                Lk = L_vector[k]
                # im - min jn - medium kl - max
                index1 = ((i -1) * M + (j ) -1)*M + (k) # i<j<k
                index2 = ((i -1) * M + (k ) -1)*M + (j) # i<k<j
                index3 = ((j -1) * M + (i ) -1)*M + (k) # j<i<K
                index4 = ((k -1) * M + (i ) -1)*M + (j) # j<k<i
                index5 = ((j -1) * M + (k ) -1)*M + (i) # k<i<j
                index6 = ((k -1) * M + (j ) -1)*M + (i) # k<j<i
                coeff_matrix_phi_thirdOrder[index1, index1] = -1*(1/Li + 1/Lj + 1/Lk)
                coeff_matrix_phi_thirdOrder[index2, index2] = 1
                coeff_matrix_phi_thirdOrder[index3, index3] = 1
                coeff_matrix_phi_thirdOrder[index4, index4] = 1
                coeff_matrix_phi_thirdOrder[index5, index5] = 1
                coeff_matrix_phi_thirdOrder[index6, index6] = 1
                b_phi_thirdOrder[index1] = -2*mutation_rate/(1+mutation_rate) * ( phi_secondOrder_matrix[ j, k ]/(numOfStrategies * Li) + phi_secondOrder_matrix[ i, k ]/(numOfStrategies * Lj) + phi_secondOrder_matrix[i, j ]/(numOfStrategies * Lk))
                b_phi_thirdOrder[index2] = 0
                b_phi_thirdOrder[index3] = 0
                b_phi_thirdOrder[index4] = 0
                b_phi_thirdOrder[index5] = 0
                b_phi_thirdOrder[index6] = 0
                coeff_matrix_phi_thirdOrder[index2, index1] = -1 
                coeff_matrix_phi_thirdOrder[index3, index1] = -1 
                coeff_matrix_phi_thirdOrder[index4, index1] = -1 
                coeff_matrix_phi_thirdOrder[index5, index1] = -1 
                coeff_matrix_phi_thirdOrder[index6, index1] = -1 

                for h in 1:M
                    coeff_matrix_phi_thirdOrder[index1, ((h -1) * M + (j ) -1)*M + (k)] += (1-mutation_rate)/(1+mutation_rate) * A[i, h] / Li 
                    coeff_matrix_phi_thirdOrder[index1, ((i -1) * M + (h ) -1)*M + (k)] += (1-mutation_rate)/(1+mutation_rate) * A[j, h] / Lj 
                    coeff_matrix_phi_thirdOrder[index1, ((i -1) * M + (j ) -1)*M + (h)] += (1-mutation_rate)/(1+mutation_rate) * A[k, h] / Lk 
                end
            end
        end
    end

    for i in 1:M
        for j in 1:M
            index1 = ((i -1) * M + (i ) -1)*M + (j) # i=j!=k
            index2 = ((i -1) * M + (j ) -1)*M + (i) # i=k!=j
            index3 = ((i -1) * M + (j ) -1)*M + (j) # k=j!=i
            coeff_matrix_phi_thirdOrder[index1, index1] = 1
            coeff_matrix_phi_thirdOrder[index2, index2] = 1
            coeff_matrix_phi_thirdOrder[index3, index3] = 1
            b_phi_thirdOrder[index1] = phi_secondOrder_matrix[ i, j ]
            b_phi_thirdOrder[index2] = phi_secondOrder_matrix[ i, j ]
            b_phi_thirdOrder[index3] = phi_secondOrder_matrix[ i, j ]
        end
        index_diag = ((i -1) * M + (i ) -1)*M + (i)
        b_phi_thirdOrder[index_diag] = 1
    end
    

    #phi_thirdOrder = coeff_matrix_phi_thirdOrder \ b_phi_thirdOrder
    phi_thirdOrder = idrs(coeff_matrix_phi_thirdOrder , b_phi_thirdOrder)
    return phi_thirdOrder
end

function multiple_calculate_m_coefficient_PC(A, I, num_nodes, offset, group_list, alpha, k,l,i,m,j,n)
    # default: alpha==0
    if (i==j)&&(m==n)
        sum_A_jn_kg = 0
        for g in 1:group_list[k]
            sum_A_jn_kg += A[n+offset[j], g+offset[k]]
        end
        m_kl_im_jn=(I[j,k] - sum_A_jn_kg)/(4*num_nodes*group_list[j])

    else
        m_kl_im_jn= A[n+offset[j], m+offset[i]]*(I[i,k] - I[j,k])/(4*num_nodes*group_list[j])
    end

    return m_kl_im_jn
end

function multiple_calculate_m_coefficient_PC_vector(A, I, num_nodes, offset, group_list)
    M = sum(group_list)
    m_klimjn = zeros(M^3)
    for i in 1:num_nodes
        for m in 1:group_list[i]
            for j in 1:num_nodes
                for n in 1:group_list[j]
                    for k in 1:num_nodes
                        for l in 1:group_list[k]
                            # default: alpha==0
                            if (i==j)&&(m==n)
                                sum_A_jn_kg = 0
                                for g in 1:group_list[k]
                                    sum_A_jn_kg += A[n+offset[j], g+offset[k]]
                                end
                                m_klimjn[ ((l + offset[k] -1) * M + (m + offset[i] ) -1)*M + (n + offset[j]) ]=(I[j,k] - sum_A_jn_kg)/(4*num_nodes*group_list[j])

                            else
                                m_klimjn[ ((l + offset[k] -1) * M + (m + offset[i] ) -1)*M + (n + offset[j]) ]= A[n+offset[j], m+offset[i]]*(I[i,k] - I[j,k])/(4*num_nodes*group_list[j])
                            end
                        end
                    end
                end
            end
        end
    end

    return m_klimjn
end

function multiple_donationGame_K1K2_PC(edge_list, num_nodes, group_list, group_data, mutation_rate, alpha=0)
    offset = calculate_offset(group_list)
    D = calculate_D(edge_list, num_nodes, group_list, group_data, offset)
    A = calculate_A(D)
    adj = edgeList2adjacencyMatrix(edge_list, num_nodes)
    W = adj2highorder_adj(adj, num_nodes, group_list, offset, group_data)
    M = sum(group_list)
    omega_matrix = multiple_calculate_omega(W, num_nodes, offset, group_list)
    I = Diagonal(ones(num_nodes)) 

    K1 = 0
    K2 = 0

    if mutation_rate > 0

        pi_mut = multiple_calculate_pi_mut_PC(A, num_nodes, offset, group_list, mutation_rate)
        #@time begin
        phi_secondOrder_matrix = multiple_calculate_phi_secondOrder_PC(A, num_nodes, offset, group_list, mutation_rate, 2)
        #end
        #println("phi second order:$(phi_secondOrder_matrix)")
        #@time begin
        for i in 1:num_nodes
            for m in 1:group_list[i]
                for j in 1:num_nodes
                    for n in 1:group_list[j]
                        for k in 1:num_nodes
                            for l in 1:group_list[k]
                                for h in 1:num_nodes
                                    for s in 1:group_list[h]
                                        K1 += 1/mutation_rate * pi_mut[m + offset[i]] * multiple_calculate_m_coefficient_PC(A, I, num_nodes, offset, group_list, alpha, k,l,j,n,i,m) * omega_matrix[l + offset[k], s + offset[h]] * ( 2*(1-mutation_rate)*phi_secondOrder_matrix[n+offset[j], l + offset[k]] - 2 *phi_secondOrder_matrix[m + offset[i], l + offset[k]] + mutation_rate )
                                        K2 += 1/mutation_rate * pi_mut[m + offset[i]] * multiple_calculate_m_coefficient_PC(A, I, num_nodes, offset, group_list, alpha, k,l,j,n,i,m) * omega_matrix[l + offset[k], s + offset[h]] * ( 2*(1-mutation_rate)*phi_secondOrder_matrix[n+offset[j], s + offset[h]] - 2 *phi_secondOrder_matrix[m + offset[i], s + offset[h]] + mutation_rate )
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
        #end     
    else
        println("Wrong mutation rate")
        return 0
    end

    return K1, K2

end


function multiple_calculate_frequencyOfStrategy_vector_structureCoefficientKonwn(PayoffMatrix, numOfStrategies, lambda1, lambda2, lambda3, selection_intensity)
    a_yStar_bar, a_starY_bar, a_starStar_bar, a_bar = multiple_payoffProcessing(PayoffMatrix, numOfStrategies)
    
    frequency_vector = zeros(numOfStrategies)

    for y in 1:numOfStrategies
        frequency_vector[y] = 1/numOfStrategies + selection_intensity * 1/numOfStrategies * (lambda1 * (PayoffMatrix[y,y] - a_starStar_bar) +lambda2 * (a_yStar_bar[y] - a_starY_bar[y]) + lambda3*( a_yStar_bar[y] - a_bar))
    end

    return frequency_vector
end