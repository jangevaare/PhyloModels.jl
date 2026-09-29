# default_rng was introduced after Julia 1.0; keep the supported old runtime.
@static if isdefined(Random, :default_rng)
  _default_rng() = Random.default_rng()
else
  _default_rng() = Random.GLOBAL_RNG
end

"""
    rand([rng], DNASeq, model, n)
    rand([rng], RNASeq, model, n)

Generate a sequence from the model's stationary frequencies. Supply an
`AbstractRNG`, for example `MersenneTwister(123)`, to control the random stream.
Omitting it uses Julia's default RNG. Exact streams may differ across Julia and
dependency versions.
"""
rand(::Type{T}, mod::NASM, n::Int64) where {T <: Union{DNASeq, RNASeq}} =
  rand(_default_rng(), T, mod, n)

"""
    simulate([rng], sequence_type, tree, model, n)
    simulate([rng], sequence_type, tree, model, site_rates)

Generate sequences at every node, using the same RNG for the root and all
branches. For example, `simulate(MersenneTwister(123), DNASeq, tree, model, 100)`.
Calls without an RNG use Julia's default RNG.
"""
simulate(::Type{T}, tree::Tree, mod::NASM, n::Int64) where {T <: Union{DNASeq, RNASeq}} =
  simulate(_default_rng(), T, tree, mod, n)
simulate(::Type{T}, tree::Tree, mod::NASM, rates::Vector{Float64}) where {T <: Union{DNASeq, RNASeq}} =
  simulate(_default_rng(), T, tree, mod, rates)

"""
    simulate!([rng], root_sequence, tree, model, site_rates)

Simulate descendants from a supplied root sequence using the specified RNG, or
Julia's default RNG when omitted. The returned dictionary retains the supplied
root sequence; the root sequence itself is not mutated.
"""
simulate!(root::Union{DNASeq, RNASeq}, tree::Tree, mod::NASM, rates::Vector{Float64}) =
  simulate!(_default_rng(), root, tree, mod, rates)

# GeneticBitArrays' weighted rand methods do not accept an RNG. Keep the
# RNG-aware sampling local rather than extending methods on foreign types.
function _rand_sequence(rng::AbstractRNG, ::Type{T}, weights::Weights, n::Int) where {T <: Union{DNASeq, RNASeq}}
  data = falses(4, n)
  for j in 1:n
    data[sample(rng, 1:4, weights), j] = true
  end
  return T(data; checkinput=false)
end

function _rand_sequence(rng::AbstractRNG, ::Type{T}, weights::AbstractVector{<:Weights}) where {T <: Union{DNASeq, RNASeq}}
  data = falses(4, length(weights))
  for j in eachindex(weights)
    data[sample(rng, 1:4, weights[j]), j] = true
  end
  return T(data; checkinput=false)
end

function rand(rng::AbstractRNG, ::Type{T}, mod::S, n::Int64) where {T<:Union{DNASeq, RNASeq}, S <: NASM}
  return _rand_sequence(rng, T, Weights(_π(mod)), n)
end

function simulate!(rng::AbstractRNG, root_seq::T,
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
    node_data[i] = _rand_sequence(rng, T, wv)
  end
  return node_data
end


function simulate(rng::AbstractRNG, ::Type{T},
                  tree::Tree,
                  mod::S,
                  site_rates::Vector{Float64}) where {T <: Union{DNASeq, RNASeq}, S <: NASM}
  # Simulation order
  visit_order = reverse(postorder(tree))
  # Sequence length
  len = length(site_rates)
  node_data = Dict{Int64, T}()
  # Generate root sequence
  node_data[visit_order[1]] = rand(rng, T, mod, len)
  # Iterate through nodes
  for i in visit_order[2:end]
    source = tree.branches[tree.nodes[i].in[1]].source
    source_seq = node_data[source]
    branch_length = tree.branches[tree.nodes[i].in[1]].length
    wv = _transition_weights(mod, branch_length, source_seq.data, site_rates)
    node_data[i] = _rand_sequence(rng, T, wv)
  end
  return node_data
end


function simulate(rng::AbstractRNG, ::Type{T},
                   tree::Tree,
                   mod::S,
                   n::Int64) where {T <: Union{DNASeq, RNASeq}, S <: NASM}
  # Simulation order
  visit_order = reverse(postorder(tree))
  node_data = Dict{Int64, T}()
  # Generate root sequence
  node_data[visit_order[1]] = rand(rng, T, mod, n)
  # Iterate through nodes
  for i in visit_order[2:end]
    source = tree.branches[tree.nodes[i].in[1]].source
    source_seq = node_data[source]
    branch_length = tree.branches[tree.nodes[i].in[1]].length
    pmat = P(mod, branch_length)
    wv = _transition_weights(pmat, source_seq.data)
    node_data[i] = _rand_sequence(rng, T, wv)
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
