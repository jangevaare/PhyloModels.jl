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
