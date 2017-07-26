# PhyloModels.jl

[![Build Status](https://travis-ci.org/jangevaare/PhyloModels.jl.svg?branch=master)](https://travis-ci.org/jangevaare/PhyloModels.jl)

PhyloModels.jl is a package for performing phylogenetic simulation and inference in Julia. This package is no longer being developed or maintained in favour of [BioJulia/SubstitutionModels.jl](https://github.com/BioJulia/SubstitutionModels.jl) and related bioinformatics packages in the [BioJulia](https://github.com/BioJulia) Ecosystem.


### Simulation
Genetic sequences can be simulated from phylogenetic trees by defining a root sequence, site rates, and a nucleotide substitution model. The following nucleotide substitution models are currently available: `JC69`, `K80`, `F81`, `F84`, `HKY85`, `TN93`, and `GTR`.


### Inference
The log likelihood of phylogenetic trees can be calculated when genetic sequences have been observed at all leaves, and a nucleotide substitution model has been specified by using the `loglikelihood` function. Below is an example from *Molecular Evolution: A Statistical Approach*

    > using PhyloTrees, PhyloModels
    >
    > # Build tree
    > tree = Tree()
    > addnodes!(tree, 9)
    > addbranch!(tree, 9, 6, 0.1)
    > addbranch!(tree, 9, 8, 0.1)
    > addbranch!(tree, 6, 7, 0.1)
    > addbranch!(tree, 6, 3, 0.2)
    > addbranch!(tree, 7, 1, 0.2)
    > addbranch!(tree, 7, 2, 0.2)
    > addbranch!(tree, 8, 4, 0.2)
    > addbranch!(tree, 8, 5, 0.2)

    Phylogenetic tree with 9 nodes and 8 branches

    > # Set state of leaf nodes
    > node_data = Dict{Int64, Sequence}()
    > node_data[1] = Sequence("T")
    > node_data[2] = Sequence("C")
    > node_data[3] = Sequence("A")
    > node_data[4] = Sequence("C")
    > node_data[5] = Sequence("C")
    >
    > # Parametrize substitution model
    > model = K80([2.])

    Kimura 1980 substitution model

    > # Calculate log likelihood
    > loglikelihood(tree, model, node_data)

    -7.5814075725577
