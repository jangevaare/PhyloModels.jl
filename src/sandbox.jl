type PhyloTrace
  substitutionmodel::Vector{SubstitutionModel}
  tree::Vector{Tree}
  logposterior::Vector{Float64}
end


type PhyloIteration
  substitutionmodel::SubstitutionModel
  tree::Tree
  logposterior::Float64
end


PhyloProposal = PhyloIteration


function push!(trace::PhyloTrace, iteration::PhyloIteration)
  push!(trace.substitutionmodel, iteration.substitutionmodel)
  push!(trace.tree, iteration.tree)
  push!(trace.logposterior, iteration.logposterior)
end


function append!(trace1::PhyloTrace, trace2::PhyloTrace)
  append!(trace1.substitutionmodel, trace2.substitutionmodel)
  append!(trace1.tree, trace2.tree)
  append!(trace1.logposterior, trace2.logposterior)
end


function deleteat!(trace::PhyloTrace, inds)
  deleteat!(trace.substitutionmodel, inds)
  deleteat!(trace.tree, inds)
  deleteat!(trace.logposterior, inds)
  return trace
end


function length(x::PhyloTrace)
  return length(x.substitutionmodel)
end


function show(io::IO, object::PhyloTrace)
  print(io, "PhyloTrace object (MCMC iterations: $(length(object)))")
end


function show(io::IO, object::PhyloIteration)
  print(io, "MCMC iteration object")
end


"""
propose(currentstate::SubstitutionModel,
        variance::Array{Float64})

Generate a `SubstitutionModel` proposal using the multivariate normal
distribution as the transition kernel, with a previous set of
`SubstitutionModel` parameters as the mean vector and a transition kernel
variance as the variance-covariance matrix
"""
function propose(currentstate::SubstitutionModel,
                 variance::Array{Float64})
  newstate = copy(currentstate)
  newstate.Θ = rand(MvNormal(newstate.Θ, variance))
  return newstate
end
