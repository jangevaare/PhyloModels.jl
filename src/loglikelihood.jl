"""
    LikelihoodWorkspace(tree, node_data; compress=true)

Allocate buffers for repeated likelihood evaluations of the same tree topology
and sequence length. Branch lengths, model parameters, and sequence contents may
change between calls to `loglikelihood!`. Rebuild the workspace after topology or
sequence-length changes. Identical alignment columns are evaluated once by default;
set `compress=false` to evaluate every site separately. Pattern groups are refreshed
when sequence contents change. All sites use the same model and branch lengths.
A workspace must not be shared by concurrent calls.
"""
mutable struct LikelihoodWorkspace
  tree::Tree
  order::Vector{Int64}
  leaves::Vector{Int64}
  topology::Dict{Int64, Tuple{Vector{Int64}, Vector{Int64}}}
  endpoints::Dict{Int64, Tuple{Int64, Int64}}
  calculations::Dict{Int64, Matrix{Float64}}
  scratch::Matrix{Float64}
  nsites::Int
  compress::Bool
  representatives::Vector{Int}
  site_pattern::Vector{Int}
  multiplicities::Vector{Int}
  sequences::Dict{Int64, BitMatrix}
end

function LikelihoodWorkspace(tree::Tree, node_data::Union{NodeDNA, NodeRNA}; compress::Bool=true)
  isempty(tree.nodes) && throw(ArgumentError("Tree must contain a root"))
  roots = [id for (id, node) in tree.nodes if isempty(node.in)]
  length(roots) == 1 || throw(ArgumentError("Tree must have exactly one root"))
  all(length(node.in) <= 1 for node in values(tree.nodes)) ||
    throw(ArgumentError("Nodes must have at most one parent"))
  order = postorder(tree)
  leaves = [id for id in order if isempty(tree.nodes[id].out)]
  all(haskey(node_data, id) for id in leaves) || error("Some leaves are missing sequence data")
  n = length(node_data[first(leaves)])
  all(length(node_data[id]) == n for id in leaves) ||
    throw(DimensionMismatch("Leaf sequences must have equal lengths"))
  topology = Dict(id => (copy(node.in), copy(node.out)) for (id, node) in tree.nodes)
  endpoints = Dict(id => (b.source, b.target) for (id, b) in tree.branches)
  representatives, site_pattern, multiplicities = compress ?
    _site_patterns(leaves, node_data, n) : (collect(1:n), collect(1:n), ones(Int, n))
  npatterns = length(representatives)
  calculations = Dict(id => Matrix{Float64}(undef, 4, npatterns) for id in order)
  sequences = compress ? Dict(id => copy(node_data[id].data) for id in leaves) : Dict{Int64, BitMatrix}()
  return LikelihoodWorkspace(tree, order, leaves, topology, endpoints,
                             calculations, Matrix{Float64}(undef, 4, npatterns),
                             n, compress, representatives, site_pattern,
                             multiplicities, sequences)
end

function _check_workspace(w::LikelihoodWorkspace, node_data)
  tree = w.tree
  valid = length(tree.nodes) == length(w.topology) && length(tree.branches) == length(w.endpoints)
  if valid
    for (id, (incoming, outgoing)) in w.topology
      if !haskey(tree.nodes, id) || tree.nodes[id].in != incoming || tree.nodes[id].out != outgoing
        valid = false
        break
      end
    end
    for (id, endpoints) in w.endpoints
      if !haskey(tree.branches, id) || (tree.branches[id].source, tree.branches[id].target) != endpoints
        valid = false
        break
      end
    end
  end
  valid || throw(ArgumentError("Tree topology changed; rebuild LikelihoodWorkspace"))
  for id in w.leaves
    haskey(node_data, id) || error("Some leaves are missing sequence data")
    length(node_data[id]) == w.nsites ||
      throw(DimensionMismatch("Sequence length changed; rebuild LikelihoodWorkspace"))
  end
  return nothing
end

"""
    loglikelihood!(workspace, model, node_data; output_calculations=false)

Evaluate the likelihood using reusable buffers. Sequence data are read afresh on
every call. With `output_calculations=true`, return independent copies of the
calculation matrices and visit order, as well as the log likelihood.
"""
function loglikelihood!(w::LikelihoodWorkspace, mod::NASM,
                        node_data::Union{NodeDNA, NodeRNA}; output_calculations::Bool=false)
  _check_workspace(w, node_data)
  _refresh_patterns!(w, node_data)
  for id in w.leaves
    dest = w.calculations[id]
    source = node_data[id].data
    for j in eachindex(w.representatives), k in 1:4
      dest[k, j] = source[k, w.representatives[j]]
    end
  end
  for id in w.order
    edges = w.tree.nodes[id].out
    for (k, edge) in enumerate(edges)
      branch = w.tree.branches[edge]
      p = P(mod, branch.length)
      child = w.calculations[branch.target]
      if k == 1
        mul!(w.calculations[id], p, child)
      else
        mul!(w.scratch, p, child)
        w.calculations[id] .*= w.scratch
      end
    end
  end
  root = w.calculations[last(w.order)]
  frequencies = _π(mod)
  ll = 0.0
  for j in axes(root, 2)
    probability = 0.0
    for k in 1:4
      probability += frequencies[k] * root[k, j]
    end
    ll += w.multiplicities[j] * log(probability)
  end
  if output_calculations
    return ll, Dict(id => x[:, w.site_pattern] for (id, x) in w.calculations), copy(w.order)
  end
  return ll
end

"""
    loglikelihood(tree, model, node_data; output_calculations=false, compress=true)

Compute the alignment log likelihood. Identical site patterns are compressed by
default; use `compress=false` to disable grouping. Returned calculation matrices
always contain all sites in their original order. For repeated evaluations,
construct a `LikelihoodWorkspace` once and call `loglikelihood!`.
"""
function loglikelihood(tree::Tree, mod::NASM, node_data::Union{NodeDNA, NodeRNA};
                       output_calculations::Bool=false, compress::Bool=true)
  return loglikelihood!(LikelihoodWorkspace(tree, node_data; compress=compress), mod, node_data;
                        output_calculations=output_calculations)
end


# Encode the four base flags per leaf, preserving ambiguous symbols and gaps.
# Keys are copied only for new patterns, never mutated after insertion.
function _site_patterns(leaves, node_data, n)
  patterns = Dict{Vector{UInt8}, Int}()
  key = Vector{UInt8}(undef, length(leaves))
  representatives = Int[]
  multiplicities = Int[]
  site_pattern = Vector{Int}(undef, n)
  for j in 1:n
    for (i, id) in enumerate(leaves)
      data = node_data[id].data
      mask = UInt8(0)
      for k in 1:4
        mask |= UInt8(data[k, j]) << (k-1)
      end
      key[i] = mask
    end
    pattern = get(patterns, key, 0)
    if pattern == 0
      push!(representatives, j)
      push!(multiplicities, 0)
      pattern = length(representatives)
      patterns[copy(key)] = pattern
    end
    site_pattern[j] = pattern
    multiplicities[pattern] += 1
  end
  return representatives, site_pattern, multiplicities
end

function _refresh_patterns!(w::LikelihoodWorkspace, node_data)
  w.compress || return nothing
  all(w.sequences[id] == node_data[id].data for id in w.leaves) && return nothing
  representatives, site_pattern, multiplicities = _site_patterns(w.leaves, node_data, w.nsites)
  n = length(representatives)
  if n != size(w.scratch, 2)
    w.calculations = Dict(id => Matrix{Float64}(undef, 4, n) for id in w.order)
    w.scratch = Matrix{Float64}(undef, 4, n)
  end
  w.representatives = representatives
  w.site_pattern = site_pattern
  w.multiplicities = multiplicities
  for id in w.leaves
    w.sequences[id] .= node_data[id].data
  end
  return nothing
end
