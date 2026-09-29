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
    wv = _transition_weights(mod, branch_length, source_seq.data, site_rates)
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
    wv = _transition_weights(mod, branch_length, source_seq.data, site_rates)
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
    wv = _transition_weights(pmat, source_seq.data)
    node_data[i] = rand(T, wv, checkinput=false)
  end
  return node_data
end


# Sampling only reads these weights, so sites with the same unambiguous source
# nucleotide can share a distribution. Ambiguities retain the matrix product.
_weight_table(p) = [Weights(collect(p[:, k])) for k in 1:4]

function _site_weights(p, table, data, j)
  state = 0
  for k in 1:4
    if data[k, j]
      state == 0 || return Weights(collect(p * view(data, :, j)))
      state = k
    end
  end
  return state == 0 ? Weights(collect(p * view(data, :, j))) : table[state]
end

function _transition_weights(p, data)
  table = _weight_table(p)
  return [_site_weights(p, table, data, j) for j in axes(data, 2)]
end

function _transition_weights(mod, branch_length, data, rates)
  result = Vector{typeof(Weights(zeros(4)))}(undef, length(rates))
  isempty(rates) && return result
  first_p = P(mod, branch_length * rates[1])
  counts = Dict{Float64, Int}()
  for rate in rates
    counts[rate] = get(counts, rate, 0) + 1
  end
  # Cache only repeated categories, avoiding a matrix and four distributions
  # per site for continuous rate distributions.
  cache = Dict{Float64, Tuple{typeof(first_p), Vector{eltype(result)}}}()
  for j in eachindex(rates)
    rate = rates[j]
    if counts[rate] > 1
      p, table = get!(cache, rate) do
        p = j == 1 ? first_p : P(mod, branch_length * rate)
        (p, _weight_table(p))
      end
      result[j] = _site_weights(p, table, data, j)
    else
      p = j == 1 ? first_p : P(mod, branch_length * rate)
      result[j] = Weights(collect(p * view(data, :, j)))
    end
  end
  return result
end
