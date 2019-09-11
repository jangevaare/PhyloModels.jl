module PhyloModels

  using PhyloTrees,
        SubstitutionModels,
        GeneticBitArrays

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

  include("core.jl")

  export NodeDNA, loglikelihood

end # module
