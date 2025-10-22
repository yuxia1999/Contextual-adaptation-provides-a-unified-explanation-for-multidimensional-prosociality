using LinearAlgebra
using SparseArrays
using Printf
using Random
using Distributions
using StatsBase

function split_edges(edge_list::Array{Tuple{Int, Int, Int}, 1}, num_nodes::Int, M::Int, method::Symbol, gamma::Float64=1.0)
    out_degree = Dict{Int, Int}()
    for i in 0:num_nodes-1
        out_degree[i] = 0
    end
    
    for edge in edge_list
        source = edge[1]
        out_degree[source] += 1
    end
    
    L = Dict{Int, Int}()
    remaining_groups = M
    
    for i in 0:num_nodes-1
        L[i] = 1
        remaining_groups -= 1
    end
    
    if method == :weighted_random
        weights = [out_degree[i]-1 for i in 0:num_nodes-1] 
        total_weight = sum(weights)
        
        while remaining_groups > 0
            probs = [weights[i] / total_weight for i in 1:num_nodes]  
            selected_node = rand(Categorical(probs)) - 1  
            if L[selected_node] < out_degree[selected_node]
                L[selected_node] += 1
                remaining_groups -= 1
                weights[selected_node + 1] -= 1  
                total_weight -= 1
            end
        end
        
    elseif method == :positive_infinite
        sorted_nodes = sort(collect(0:num_nodes-1), by=i -> out_degree[i], rev=true)  
        for i in sorted_nodes
            while remaining_groups > 0 && L[i] < out_degree[i]
                L[i] += 1
                remaining_groups -= 1
            end
        end
    
    elseif method == :negative_infinite
        sorted_nodes = sort(collect(0:num_nodes-1), by=i -> out_degree[i])  
        for i in sorted_nodes
            while remaining_groups > 0 && L[i] < out_degree[i]
                L[i] += 1
                remaining_groups -= 1
            end
        end

    elseif method == :weighted_deterministic
        weights = [(out_degree[i]-1)^gamma for i in 0:num_nodes-1]
        total_weight = sum(weights)
        
        expected_groups = [remaining_groups * weights[i+1]/ total_weight for i in 0:num_nodes-1]
        List_done = Set{Int}()

        while true
            updated = false
            for i in 0:num_nodes-1
                if i in List_done
                    continue
                end
    
                if expected_groups[i+1] > out_degree[i] - 1
                    expected_groups[i+1] = out_degree[i] - 1
                    push!(List_done, i)
                    updated = true
                end
            end
    
            if !updated
                break
            end
    
            expected_remaining = remaining_groups - sum(expected_groups[i+1] for i in List_done)
            for i in List_done
                weights[i+1] = 0
            end
    
            total_weight = sum(weights)
            for i in 0:num_nodes-1
                if i in List_done
                    continue
                end
                expected_groups[i+1] = expected_remaining * weights[i+1] / total_weight
            end
        end

        for i in 0:num_nodes-1
            allocated_groups = floor(Int, expected_groups[i+1])
            L[i] += allocated_groups 
            remaining_groups -= allocated_groups
        end

        weights_new = [(expected_groups[i+1]-L[i]+1) for i in 0:num_nodes-1]
        total_weight_new = sum(weights_new)
        if total_weight_new !=0
            probs = weights_new / total_weight_new  

            while remaining_groups > 0
                selected_node = rand(Categorical(probs)) - 1  
                if L[selected_node] < (expected_groups[selected_node+1]+1)
                    L[selected_node] += 1
                    remaining_groups -= 1
                end
            end
        end

    elseif method == :weighted_absolute_deterministic
        weights = [(out_degree[i]-1)^gamma for i in 0:num_nodes-1]
        total_weight = sum(weights)
        
        expected_groups = [remaining_groups * weights[i+1]/ total_weight for i in 0:num_nodes-1]
        List_done = Set{Int}()

        while true
            updated = false
            for i in 0:num_nodes-1
                if i in List_done
                    continue
                end
    
                if expected_groups[i+1] > out_degree[i] - 1
                    expected_groups[i+1] = out_degree[i] - 1
                    push!(List_done, i)
                    updated = true
                end
            end
    
            if !updated
                break
            end
    
            expected_remaining = remaining_groups - sum(expected_groups[i+1] for i in List_done)
            for i in List_done
                weights[i+1] = 0
            end
    
            total_weight = sum(weights)
            for i in 0:num_nodes-1
                if i in List_done
                    continue
                end
                expected_groups[i+1] = expected_remaining * weights[i+1] / total_weight
            end
        end
    
        for i in 0:num_nodes-1
            allocated_groups = floor(Int, expected_groups[i+1])
            L[i] += allocated_groups 
            remaining_groups -= allocated_groups
        end

        weights_new = [(expected_groups[i+1]-L[i]+1) for i in 0:num_nodes-1]
        total_weight_new = sum(weights_new)
        
        if total_weight_new != 0
            probs = weights_new / total_weight_new  

            sorted_indices = sortperm(probs, rev=true)  
        
            for idx in sorted_indices
                if L[idx-1] < (expected_groups[idx] + 1) && remaining_groups > 0
                    L[idx-1] += 1
                    remaining_groups -= 1
                end
                if remaining_groups == 0
                    break
                end
            end
        end

    elseif method == :uniform
        nodes = shuffle(collect(0:num_nodes-1))
        while remaining_groups > 0
            for i in nodes
                if remaining_groups == 0
                    break
                end
                if L[i] < out_degree[i]
                    L[i] += 1
                    remaining_groups -= 1
                end
            end
        end
    else
        error("Unknown method: $method. Use :uniform or :weighted.")
    end
    
    edge_groups = Dict{Int, Vector{Vector{Int}}}()
    for i in 0:num_nodes-1
        edge_groups[i] = [Vector{Int}() for _ in 1:L[i]]
    end
    
    for i in 0:num_nodes-1
        out_edges = shuffle([edge[2] for edge in edge_list if edge[1] == i])  
        num_groups = L[i]
        group_size = length(out_edges) ÷ num_groups
        remainder = length(out_edges) % num_groups
        
        start = 1
        for j in 1:num_groups
            size = group_size + (j <= remainder ? 1 : 0)
            edge_groups[i][j] = out_edges[start:start+size-1]
            start += size
        end
    end
    
    group_list = [L[i] for i in 0:num_nodes-1]
    group_data = [edge_groups[i] for i in 0:num_nodes-1]
    
    return group_list, group_data
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

# network generation

function generate_and_export_power_graph(n::Int, k::Int, v::Real)
    if k < n-1
        @assert k % 2 == 0 "k must be even"
        @assert n >= k+1    "n must be ≥ k+1"

        adj = [Int[] for _ in 1:n]       
        deg = zeros(Int, n)              

        m = k + 1
        for i in 1:m-1, j in i+1:m
            push!(adj[i], j); push!(adj[j], i)
            deg[i] += 1; deg[j] += 1
        end

        halfk = k ÷ 2
        for new_node in (m+1):n
            existing = collect(1:new_node-1)
            weights = [deg[u]^v for u in existing]  

            targets = Int[]
            idxs = copy(existing)
            ws   = copy(weights)
            for _ in 1:halfk
                cum = cumsum(ws)
                r = rand() * cum[end]
                sel = searchsortedfirst(cum, r)
                push!(targets, idxs[sel])
                deleteat!(idxs, sel)
                deleteat!(ws,   sel)

            end

            for t in targets
                push!(adj[new_node], t)
                push!(adj[t],        new_node)
                deg[new_node] += 1
                deg[t]        += 1
            end
        end
    elseif k == n-1
        adj = [Int[] for _ in 1:n]       
        deg = zeros(Int, n)              
        m = n
        for i in 1:m-1, j in i+1:m
            push!(adj[i], j); push!(adj[j], i)
            deg[i] += 1; deg[j] += 1
        end
    else
        error("k must be less than n-1")
    end

    edge_list = Vector{Tuple{Int, Int, Int}}()
    for u in 1:n, v2 in adj[u]
        push!(edge_list, (u-1, v2-1, 1))
    end

    return edge_list

end