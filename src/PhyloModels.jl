__precompile__()

module PhyloModels

# Dependencies
using
  PhyloTrees,
  Distributions

# Methods expanded
import
  Base.show,
  Base.push!,
  Base.append!,
  Base.length,
  Base.getindex,
  Base.rand,
  Base.copy,
  Base.deleteat!,
  StatsBase.loglikelihood

export
  # Sequences
  Sequence,

  # Substitution Models
  SubstitutionModel,
  SubstitutionModelPrior,
  JC69,
  K80,
  F81,
  F84,
  HKY85,
  TN93,
  GTR,
  JC69Prior,
  K80Prior,
  F81Prior,
  F84Prior,
  HKY85Prior,
  TN93Prior,
  GTRPrior,
  Q,
  P,
  logprior,

  # Simulation
  simulate,
  simulate!,

  # Loglikelihoods
  loglikelihood

  # Package files
  include("sequences.jl")
  include("substitution_models/abstract.jl")
  include("substitution_models/jc69.jl")
  include("substitution_models/k80.jl")
  include("substitution_models/f81.jl")
  include("substitution_models/f84.jl")
  include("substitution_models/hky85.jl")
  include("substitution_models/tn93.jl")
  include("substitution_models/gtr.jl")
  include("simulation.jl")
  include("loglikelihoods.jl")

end # module
