"""
simulate(n::Int64, mod::SubstitutionModel)

Simulate a sequence of length `n`, using the base frequencies of a specified
substitution model
"""
function simulate(n::Int64,
                  mod::SubstitutionModel)
  return Sequence(convert(Array{Bool, 2}, rand(Multinomial(1, mod.π), n)))
end



"""
simulate!{B<:Any}(tree::Tree{Sequence, B},
                  mod::SubstitutionModel,
                  site_rates::Vector{Float64})

Simulate sequences for all nodes in a phylogenetic `tree` following a specified
substitution model
"""
function simulate!{B<:Any}(tree::Tree{Sequence, B},
                           mod::SubstitutionModel,
                           site_rates::Vector{Float64})
  # Simulation order
  visit_order = reverse(postorder(tree))

  # Error checking
  if !hasdata(tree.nodes[visit_order[1]])
    error("Root sequence must be specificed (node $(visit_order[1]))")
  elseif length(getdata(tree.nodes[visit_order[1]])) != length(site_rates)
    error("Dimension of root sequence must match length of site rates")
  end

  # Iterate through remaining nodes
  for i in visit_order[2:end]
    source = tree.branches[tree.nodes[i].in[1]].source
    source_seq = getdata(tree.nodes[source])
    branch_length = tree.branches[tree.nodes[i].in[1]].length
    seq = [findfirst(rand(Multinomial(1, (source_seq.nucleotides[:, j]' * P(mod, branch_length * site_rates[j]))[:]))) for j = 1:length(source_seq)]
    setdata!(tree.nodes[i], Sequence(seq))
  end
  return tree
end
