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


@testset "Identical alignment columns" begin
  for seqtype in (DNASeq, RNASeq)
    data = Dict(id => seqtype("ACNACN") for id in findleaves(tree))
    w = LikelihoodWorkspace(tree, data)
    @test size(w.scratch, 2) == 3
    @test w.multiplicities == [2, 2, 2]
    @test loglikelihood!(w, K80(2.0), data) ≈ reference_loglikelihood(tree, K80(2.0), data)
    ll, matrices, order = loglikelihood!(w, K80(2.0), data; output_calculations=true)
    oldll, oldmatrices, _ = reference_loglikelihood(tree, K80(2.0), data; output_calculations=true)
    @test ll ≈ oldll
    @test all(matrices[id] ≈ oldmatrices[id] for id in order)
    @test loglikelihood(tree, K80(2.0), data; compress=false) ≈ ll
    # Change one occurrence of a repeated column, in place, and split its group.
    id = first(findleaves(tree))
    data[id].data[:, 4] .= [false, true, false, false]
    @test loglikelihood!(w, JC69(), data) ≈ reference_loglikelihood(tree, JC69(), data)
    @test size(w.scratch, 2) == 4
    # Replacement sequences can merge all patterns again.
    for id in keys(data)
      data[id] = seqtype("AAAAAA")
    end
    @test loglikelihood!(w, JC69(), data) ≈ reference_loglikelihood(tree, JC69(), data)
    @test size(w.scratch, 2) == 1
    full = LikelihoodWorkspace(tree, data; compress=false)
    @test size(full.scratch, 2) == 6
  end
end


using Random
@testset "Explicit random number generators" begin
  generators = [MersenneTwister]
  if isdefined(Random, :Xoshiro)
    push!(generators, Random.Xoshiro)
  end
  rates = repeat([0.0, 0.5, 1.0, 1.0], 16)
  for generator in generators, seqtype in (DNASeq, RNASeq)
    root = seqtype(repeat("ACNR", 16))
    saved = copy(root.data)
    calls = (
      rng -> rand(rng, seqtype, K80(2.0), 64),
      rng -> simulate(rng, seqtype, tree, K80(2.0), 64),
      rng -> simulate(rng, seqtype, tree, K80(2.0), rates),
      rng -> simulate!(rng, root, tree, K80(2.0), rates),
    )
    for call in calls
      @test call(generator(123)) == call(generator(123))
      rng = generator(123)
      call(rng)
      @test rand(rng) != rand(generator(123))
      # An explicit RNG must not consume or reseed the default random stream.
      Random.seed!(987)
      expected = rand()
      Random.seed!(987)
      call(generator(123))
      @test rand() == expected
    end
    @test root.data == saved
    result = simulate!(generator(123), root, tree, K80(2.0), rates)
    @test result[last(postorder(tree))] === root
    @test isempty(rand(generator(123), seqtype, JC69(), 0).data)
    @test length(simulate(generator(123), seqtype, tree, JC69(), Float64[])[1]) == 0
    @test_throws ErrorException simulate!(generator(123), root, tree, JC69(), [1.0])
    if generator == MersenneTwister
      for call in (
        () -> rand(seqtype, K80(2.0), 64),
        () -> simulate(seqtype, tree, K80(2.0), 64),
        () -> simulate(seqtype, tree, K80(2.0), rates),
        () -> simulate!(root, tree, K80(2.0), rates),
      )
        Random.seed!(123)
        first_result = call()
        Random.seed!(123)
        @test call() == first_result
      end
    end
  end
end
