module PhyloModels

  using PhyloTrees,
        SubstitutionModels,
        GeneticBitArrays

  import SubstitutionModels._π
  import Base.rand

  # Re-export all of PhyloTrees, SubstitutionModels, and GeneticBitArrays
  for name in names(PhyloTrees)
    @eval export $(name)
  end

  for name in names(SubstitutionModels)
    @eval export $(name)
  end

  for name in names(GeneticBitArrays)
    @eval export $(name)
  end

  const NodeDNA = Dict{Int64, DNASeq}
  const NodeRNA = Dict{Int64, RNASeq}

  include("loglikelihood.jl")
  include("simulate.jl")

  export NodeDNA, NodeRNA, simulate!, simulate, rand, loglikelihood

end # module
