using Test,
      PhyloModels

# An example from "Molecular Evolution: A Statistical Approach" by Ziheng Yang

# Describe a Phylogenetic Tree
tree = Tree()

# Add 9 nodes
addnodes!(tree, 9)

# Connect with branches...
addbranch!(tree, 9, 6, 0.1)
addbranch!(tree, 9, 8, 0.1)
addbranch!(tree, 6, 7, 0.1)
addbranch!(tree, 6, 3, 0.2)
addbranch!(tree, 7, 1, 0.2)
addbranch!(tree, 7, 2, 0.2)
addbranch!(tree, 8, 4, 0.2)
addbranch!(tree, 8, 5, 0.2)

# Specify sequences of leaf nodes
node_data = NodeDNA()
node_data[1] = DNASeq("T")
node_data[2] = DNASeq("C")
node_data[3] = DNASeq("A")
node_data[4] = DNASeq("C")
node_data[5] = DNASeq("C")

# Specify a Nucleic Acid Substitution Model
model = K80(2.0)

# loglikelihood calculation
ll = loglikelihood(tree, model, node_data)

@test ll ≈ -7.5814075725577

# Simulation
node_data = simulate(DNASeq,
                     tree,
                     model,
                     1000)
@test length(node_data[1]) == 1000

@testset "Simulation distributions" begin
  for seqtype in (DNASeq, RNASeq)
    seq = seqtype(seqtype == DNASeq ? "ACGTNR" : "ACGUNR")
    p = P(K80(2.0), 0.2)
    old = [Weights(p * seq.data[:, j]) for j in 1:length(seq)]
    new = PhyloModels._transition_weights(p, seq.data)
    @test all(collect(old[j]) ≈ collect(new[j]) for j in eachindex(old))
    rates = [1.0, 0.0, 1.0, 2.0, 2.0, 0.5]
    old = [Weights(P(K80(2.0), 0.2*rates[j]) * seq.data[:, j]) for j in eachindex(rates)]
    new = PhyloModels._transition_weights(K80(2.0), 0.2, seq.data, rates)
    @test all(collect(old[j]) ≈ collect(new[j]) for j in eachindex(old))
    @test length(simulate(seqtype, tree, K80(2.0), rates)[1]) == length(rates)
    @test length(simulate!(seq, tree, K80(2.0), rates)[1]) == length(rates)
  end
  seq = DNASeq("AAAA")
  weights = PhyloModels._transition_weights(P(K80(2.0), 0.1), seq.data)
  @test weights[1] === weights[4]
  @test isempty(PhyloModels._transition_weights(P(K80(2.0), 0.1), DNASeq("").data))
end


function reference_loglikelihood(tree::Tree,
                       mod::T,
                       node_data::N;
                       output_calculations::Bool=false) where {T <: NASM, N <: Union{NodeDNA, NodeRNA}}
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
    return sum(log.(PhyloModels._π(mod)' * calculations[visit_order[end]])), calculations, visit_order
  else
    return sum(log.(PhyloModels._π(mod)' * calculations[visit_order[end]]))
  end
end


@testset "Likelihood workspace" begin
  for seqtype in (DNASeq, RNASeq)
    data = Dict(id => seqtype(seqtype == DNASeq ? "ACGTNRAC" : "ACGUNRAC") for id in findleaves(tree))
    w = LikelihoodWorkspace(tree, data)
    for model in (JC69(), K80(2.0), F81(0.1, 0.2, 0.3, 0.4))
      @test loglikelihood!(w, model, data) ≈ reference_loglikelihood(tree, model, data)
      @test loglikelihood(tree, model, data) ≈ reference_loglikelihood(tree, model, data)
    end
    value, matrices, order = loglikelihood!(w, K80(2.0), data; output_calculations=true)
    expected, old_matrices, old_order = reference_loglikelihood(tree, K80(2.0), data; output_calculations=true)
    @test value ≈ expected
    @test all(matrices[id] ≈ old_matrices[id] for id in order)
    saved = copy(matrices[last(order)])
    data[first(findleaves(tree))] = seqtype("AAAAAAAA")
    @test loglikelihood!(w, JC69(), data) ≈ reference_loglikelihood(tree, JC69(), data)
    @test matrices[last(order)] == saved
    data[first(findleaves(tree))] = seqtype("A")
    @test_throws DimensionMismatch loglikelihood!(w, JC69(), data)
  end
  data = Dict(id => DNASeq("") for id in findleaves(tree))
  @test loglikelihood(tree, JC69(), data) == 0.0
  w = LikelihoodWorkspace(tree, data)
  edge = first(keys(tree.branches))
  b = tree.branches[edge]
  tree.branches[edge] = PhyloTrees.Branch(b.source, b.target, 0.7)
  @test loglikelihood!(w, JC69(), data) == 0.0
  tree.branches[edge] = b
  changed = deepcopy(tree)
  w = LikelihoodWorkspace(changed, data)
  branch!(changed, first(findleaves(changed)), 0.1)
  @test_throws ArgumentError loglikelihood!(w, JC69(), data)
end
