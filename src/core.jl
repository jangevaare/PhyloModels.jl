const NodeDNA = Dict{Int64, DNABitArray}

function loglikelihood(tree::Tree,
                       mod::T,
                       node_data::NodeDNA;
                       output_calculations::Bool=false) where T <: NASM
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
      calculations[i] = node_data[i].nucleotides
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
    return sum(log.(mod.π' * calculations[visit_order[end]])), calculations, visit_order
  else
    return sum(log.(mod.π' * calculations[visit_order[end]]))
  end
end
