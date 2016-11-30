"""
Simulate sequences for all nodes in a phylogenetic tree following a specified
substitution model
"""
function simulate!{B<:Any}(tree::Tree{Sequence, B},
                           mod::SubstitutionModel,
                           site_rates::Vector{Float64})
  # Simulation order
  visit_order = reverse(postorder(tree))

  # Error checking
  if length(get(tree.nodes[visit_order[1]].data)) != length(site_rates)
    error("Dimension of root sequence must match length of site rates")
  end

  # Iterate through remaining nodes
  for i in visit_order[2:end]
    source = tree.branches[tree.nodes[i].in[1]].source
    branch_length = get(tree.branches[tree.nodes[i].in[1]].length)
    tree.nodes[i].data = tree.nodes[source].data
    for j in 1:length(get(tree.nodes[i].data))
      site_rate = site_rates[j]
      p = P(mod, branch_length * site_rate)
      tree.nodes[i].data.value.nucleotides[:,j] = rand(Multinomial(1, (tree.nodes[source].data.value.nucleotides[:,j]' * p)[:]))
    end
  end
  return tree
end
