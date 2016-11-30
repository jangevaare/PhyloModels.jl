type PhyloTrace
  substitutionmodel::Vector{SubstitutionModel}
  tree::Vector{Tree}
  logposterior::Vector{Float64}
  acceptance::Vector{Bool}
end


type PhyloIteration
  substitutionmodel::SubstitutionModel
  tree::Tree
  logposterior::Float64
  acceptance::Bool
end


PhyloProposal = PhyloIteration


function push!(trace::PhyloTrace, iteration::PhyloIteration)
  push!(trace.substitutionmodel, iteration.substitutionmodel)
  push!(trace.tree, iteration.tree)
  push!(trace.logposterior, iteration.logposterior)
  push!(trace.acceptance, iteration.acceptance)
end


function append!(trace1::PhyloTrace, trace2::PhyloTrace)
  append!(trace1.substitutionmodel, trace2.substitutionmodel)
  append!(trace1.tree, trace2.tree)
  append!(trace1.logposterior, trace2.logposterior)
  append!(trace1.acceptance, trace2.acceptance)
end


function length(x::PhyloTrace)
  return length(x.substitutionmodel)
end


function show(io::IO, object::PhyloTrace)
  print(io, "PhyloTrace object (MCMC iterations: $(length(object)), acceptance rate: $(trunc(sum(object.acceptance)*100/length(object), 4))%)")
end


function show(io::IO, object::PhyloIteration)
  print(io, "MCMC iteration object")
end


"""
Generate the variance-covariance matrix for a MvNormal transition kernel based
upon prior distributions
"""
function transition_kernel_variance(x::SubstitutionModelPrior)
  diagonal = Float64[]
  for i in x.Θ
    push!(diagonal, var(i)*2.38^2)
  end
  return diagonal
end


"""
Adapt the variance-covariance matrix for a MvNormal transition kernel for
`SubstitutionModel`
"""
function transition_kernel_variance(x::Vector{SubstitutionModel})
  covariance_matrix = cov([x[i].Θ for i = 1:length(x)])
  kernel_var = covariance_matrix * (2.38^2) / size(covariance_matrix, 1)
  if size(kernel_var, 1) > 1
    return diag(kernel_var)
  else
    return kernel_var
  end
end


"""
Generate a `SubstitutionModel` proposal using the multivariate normal
distribution as the transition kernel, with a previous set of
`SubstitutionModel` parameters as the mean vector and a transition kernel
variance as the variance-covariance matrix
"""
function propose(currentstate::SubstitutionModel,
                 substitutionmodel_prior::SubstitutionModelPrior,
                 variance::Vector{Float64})
  newstate = currentstate
  for i in 1:length(substitutionmodel_prior.Θ)
    lb = support(substitutionmodel_prior.Θ[i]).lb
    ub = support(substitutionmodel_prior.Θ[i]).ub
    newstate.Θ[i] = rand(Truncated(Normal(currentstate.Θ[i], variance[i]), lb, ub))
  end
  if length(fieldnames(substitutionmodel_prior)) == 2
    newstate.π = rand(Dirichlet([5; 5; 5; 5]))
  end
  return newstate
end


"""
Subtree-Prune-Regraft operator
"""
function spr(tree::Tree)
  subtreeroot = sample(findnonroots(tree))
  subtreenodes = [subtreeroot; descendantnodes(tree, subtreeroot)]
  sampleorder = sample(1:length(tree.nodes), length(tree.nodes))
  for i in sampleorder
    if !(i in subtreenodes)
      reattachmentnode = i
      break
    end
  end
  return changesource!(tree, tree.nodes[subtreeroot].in, reattachmentnode)
end


"""
Calculates the log likelihood of a tree with sequences observed at all leaves
* does not assume equal site rates
* not yet ready for production use
"""
function loglikelihood(seq::Vector{Sequence},
                       tree::Tree,
                       mod::SubstitutionModel,
                       site_rates::Vector{Float64})
  seq_length = length(site_rates)
  leaves = findleaves(tree)
  if length(leaves) !== length(seq)
    error("Number of leaves and number of observed sequences do not match")
  end
  visit_order = postorder(tree)
  seq_array = fill(1., (4, seq_length, length(tree.nodes)))
  leafindex = 0
  for i in visit_order
    if isleaf(tree.nodes[i])
      leafindex += 1
      seq_array[:, :, i] = seq_array[:, :, i] .* seq[leafindex].nucleotides
    else
      branches = tree.nodes[i].out
      for j in branches
        branch_length = get(tree.branches[j].length)
        child_node = tree.branches[j].target
        @simd for k in 1:seq_length
          seq_array[:, k, i] = seq_array[:, k, i] .* (seq_array[:, k, child_node]' * P(mod, branch_length * site_rates[k]))[:]
        end
      end
    end
  end
  ll = 0.
  @simd for i in 1:seq_length
    ll += log(sum(seq_array[:, i, visit_order[end]] .* mod.π))
  end
  return ll
end
