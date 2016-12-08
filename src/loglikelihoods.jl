"""
loglikelihood{B<:Any}(tree::Tree{Sequence, B},
                      mod::SubstitutionModel)

Calculates the log likelihood of a tree with sequences observed at all leaves
"""
function loglikelihood{B<:Any}(tree::Tree{Sequence, B},
                               mod::SubstitutionModel)
  seq_length = length(getdata(tree.nodes[findleaves(tree)[1]]))
  visit_order = postorder(tree)
  seq_array = fill(1., (4, seq_length, length(tree.nodes)))
  for i in visit_order
    if isleaf(tree.nodes[i])
      if length(getdata(tree.nodes[i])) != seq_length
        error("Unequal sequence lengths")
      end
      seq_array[:, :, i] = seq_array[:, :, i] .* getdata(tree.nodes[i]).nucleotides
    else
      branches = tree.nodes[i].out
      for j in branches
        branch_length = tree.branches[j].length
        child_node = tree.branches[j].target
        p = P(mod, branch_length)
        @simd for k in 1:seq_length
          seq_array[:, k, i] = seq_array[:, k, i] .* (seq_array[:, k, child_node]' * p)[:]
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
