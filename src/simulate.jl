function rand(::Type{T}, mod::S, n::Int64) where {T<:Union{DNASeq, RNASeq}, S <: NASM}
  return rand(T, Weights(_π(mod)), n)
end

function simulate!(root_seq::T,
                   tree::Tree,
                   mod::S,
                   site_rates::Vector{Float64}) where {T <: Union{DNASeq, RNASeq}, S <: NASM}
  # Simulation order
  visit_order = reverse(postorder(tree))
  # Sequence length
  len = length(site_rates)
  node_data = Dict{Int64, T}()
  node_data[visit_order[1]] = root_seq
  # Error checking
  if length(root_seq) != len
    throw(ErrorException("Dimension of root sequence must match length of site rates"))
  end
  # Iterate through nodes
  for i in visit_order[2:end]
    source = tree.branches[tree.nodes[i].in[1]].source
    source_seq = node_data[source]
    branch_length = tree.branches[tree.nodes[i].in[1]].length
    wv = [Weights(P(mod, branch_length * site_rates[j]) * source_seq.data[:, j]) for j = 1:len]
    node_data[i] = rand(T, wv, checkinput=false)
  end
  return node_data
end


function simulate(::Type{T},
                  tree::Tree,
                  mod::S,
                  site_rates::Vector{Float64}) where {T <: Union{DNASeq, RNASeq}, S <: NASM}
  # Simulation order
  visit_order = reverse(postorder(tree))
  # Sequence length
  len = length(site_rates)
  node_data = Dict{Int64, T}()
  # Generate root sequence
  node_data[visit_order[1]] = rand(T, mod, len)
  # Iterate through nodes
  for i in visit_order[2:end]
    source = tree.branches[tree.nodes[i].in[1]].source
    source_seq = node_data[source]
    branch_length = tree.branches[tree.nodes[i].in[1]].length
    wv = [Weights(P(mod, branch_length * site_rates[j]) * source_seq.data[:, j]) for j = 1:len]
    node_data[i] = rand(T, wv, checkinput=false)
  end
  return node_data
end


function simulate(::Type{T},
                   tree::Tree,
                   mod::S,
                   n::Int64) where {T <: Union{DNASeq, RNASeq}, S <: NASM}
  # Simulation order
  visit_order = reverse(postorder(tree))
  node_data = Dict{Int64, T}()
  # Generate root sequence
  node_data[visit_order[1]] = rand(T, mod, n)
  # Iterate through nodes
  for i in visit_order[2:end]
    source = tree.branches[tree.nodes[i].in[1]].source
    source_seq = node_data[source]
    branch_length = tree.branches[tree.nodes[i].in[1]].length
    pmat = P(mod, branch_length)
    wv = [Weights(pmat * source_seq.data[:, j]) for j = 1:n]
    node_data[i] = rand(T, wv, checkinput=false)
  end
  return node_data
end
