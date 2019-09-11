module PhyloModels

  using PhyloTrees,
        SubstitutionModels,
        GeneticBitArrays

  import SubstitutionModels._π

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
  include("loglikelihood.jl")

  export NodeDNA, loglikelihood

end # module
