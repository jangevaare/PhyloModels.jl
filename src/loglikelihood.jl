"""
    LikelihoodWorkspace(tree, node_data)

Allocate buffers for repeated likelihood evaluations of the same tree topology
and sequence length. Branch lengths, model parameters, and sequence contents may
change between calls to `loglikelihood!`. Rebuild the workspace after topology or
sequence-length changes. A workspace must not be shared by concurrent calls.
"""
struct LikelihoodWorkspace
  tree::Tree
  order::Vector{Int64}
  leaves::Vector{Int64}
  topology::Dict{Int64, Tuple{Vector{Int64}, Vector{Int64}}}
  endpoints::Dict{Int64, Tuple{Int64, Int64}}
  calculations::Dict{Int64, Matrix{Float64}}
  scratch::Matrix{Float64}
end

function LikelihoodWorkspace(tree::Tree, node_data::Union{NodeDNA, NodeRNA})
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
  calculations = Dict(id => Matrix{Float64}(undef, 4, n) for id in order)
  return LikelihoodWorkspace(tree, order, leaves, topology, endpoints,
                             calculations, Matrix{Float64}(undef, 4, n))
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
    length(node_data[id]) == size(w.scratch, 2) ||
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
  for id in w.leaves
    w.calculations[id] .= node_data[id].data
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
    ll += log(probability)
  end
  if output_calculations
    return ll, Dict(id => copy(x) for (id, x) in w.calculations), copy(w.order)
  end
  return ll
end

function loglikelihood(tree::Tree, mod::NASM, node_data::Union{NodeDNA, NodeRNA};
                       output_calculations::Bool=false)
  return loglikelihood!(LikelihoodWorkspace(tree, node_data), mod, node_data;
                        output_calculations=output_calculations)
end
