"""
loglikelihood(tree::Tree,
              mod::SubstitutionModel,
              node_data::Dict{Int64, Sequence})

Calculates the log likelihood of a tree with sequences observed at all leaves
"""
function loglikelihood(tree::Tree,
                       mod::SubstitutionModel,
                       node_data::Dict{Int64, Sequence})
  # Error checking
  if !all(map(x -> x in keys(node_data), findleaves(tree)))
    error("Some leaves are missing sequence data")
  elseif length(findroots(tree)) > 1
    error("More than one root detected")
  elseif condition

  # Create a Dict to store likelihood calculations
  calculations = Dict{Int64, Array{Float64, 2}}()
  # Find node visit order for postorder traversal
  visit_order = postorder(tree)
  for i in visit_order
    if isleaf(tree, i)
      calculations[i] = node_data[i].nucleotides
    else
      branches = tree.nodes[i].out
      for j in branches
        branch_length = tree.branches[j].length
        child_node = tree.branches[j].target
        # Calculate p matrix for specific branch length
        p = P(mod, branch_length)
        # Initialize likelihood calculation array for node if not already done
        if !haskey(calculations, i)
          calculations[i] = fill(1., size(calculations[child_node]))
        end
        # Perform likelihood calculation for each nucleotide
        @simd for k in 1:size(calculations[i], 2)
          calculations[i][:,k] .*= (calculations[child_node][:, k]' * p)[:]
        end
      end
    end
  end
  ll = 0.
  @simd for i in 1:size(calculations[visit_order[end]], 2)
    ll += log(sum(calculations[visit_order[end]][:, i] .* mod.π))
  end
  return ll
end
