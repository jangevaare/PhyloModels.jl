"""
Calculate the loglikelihood of a rooted phylogenetic tree. Genetic 
sequences for each leaf must be provided in a Dict, using Node ID as
the key
"""
function loglikelihood(tree::Tree,
                       mod::T,
                       node_data::ND;
                       output_calculations::Bool=false) where {
                       T <: NASM, 
                       N <: Union{NodeDNA, NodeRNA}, 
                       ND <: Dict{Int64, N}}
  # Error checking
  if !all(map(x -> x in keys(node_data), findleaves(tree)))
    error("Some leaves are missing sequence data")
  elseif length(findroots(tree)) > 1
    error("More than one root detected")
  end

  # Create a Dict to store likelihood calculations
  calculations = Dict{Int64, Array{Float64, 2}}()

  # Find node visit order for postorder traversal
  visit_order = postorder(tree)
  for i in visit_order
    if isleaf(tree, i)
      calculations[i] = node_data[i].data
    else
      branches = tree.nodes[i].out
      for j in branches
        branch_length = tree.branches[j].length
        child_node = tree.branches[j].target
        p = P(mod, branch_length)
        if !haskey(calculations, i)
          calculations[i] = p * calculations[child_node]
        else
          calculations[i] .*= p * calculations[child_node]
        end
      end
    end
  end
  if output_calculations
    return sum(log.(_π(mod)' * calculations[visit_order[end]])), calculations, visit_order
  else
    return sum(log.(_π(mod)' * calculations[visit_order[end]]))
  end
end
