# PhyloModels.jl

[![Build Status](https://travis-ci.org/jangevaare/PhyloModels.jl.svg?branch=master)](https://travis-ci.org/jangevaare/PhyloModels.jl)

PhyloModels.jl is a package for performing phylogenetic simulation and inference in Julia.


### Simulation
Genetic sequences can be simulated from phylogenetic trees by defining a root sequence, site rates, and a nucleotide substitution model. The following nucleotide substitution models are currently available: `JC69`, `K80`, `F81`, `F84`, `HKY85`, `TN93`, and `GTR`.


### Inference
The log likelihood of phylogenetic trees can be calculated when genetic sequences have been observed at all leaves, and a nucleotide substitution model has been specified by using the `loglikelihood` function.

### Development
PhyloModels.jl intends to be a full suite of tools for phylogenetic simulation and inference in Julia. Future versions will include
* Further support for heterogenous site rate models
* Tree operators and MCMC methods
