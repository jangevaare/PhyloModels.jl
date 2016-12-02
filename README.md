# PhyloModels.jl

[![Build Status](https://travis-ci.org/jangevaare/PhyloModels.jl.svg?branch=master)](https://travis-ci.org/jangevaare/PhyloModels.jl)

PhyloModels.jl is a package for performing phylogenetic simulation and inference in Julia.


### Simulation
Genetic sequences can be simulated from phylogenetic trees by defining a root sequence, site rates, and a nucleotide substitution model. The following nucleotide substitution models are currently available: `JC69`, `K80`, `F81`, `F84`, `HKY85`, `TN93`, and `GTR`.


### Inference
The log likelihood of phylogenetic trees can be calculated when genetic sequences have been observed at all leaves, and a nucleotide substitution model has been specified by using the `loglikelihood` function. Below is an example from *Molecular Evolution: A Statistical Approach*

    > using PhyloTrees, PhyloModels
    >
    > # Build tree
    > tree = Tree{Sequence, Void}()
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
    > leaves = findleaves(tree)
    > tree.nodes[leaves[1]].data = Nullable(Sequence("T"))
    > tree.nodes[leaves[2]].data = Nullable(Sequence("C"))
    > tree.nodes[leaves[3]].data = Nullable(Sequence("A"))
    > tree.nodes[leaves[4]].data = Nullable(Sequence("C"))
    > tree.nodes[leaves[5]].data = Nullable(Sequence("C"))
    >
    > # Parametrize substitution model
    > model = K80([2.])

    Kimura 1980 substitution model

    > # Calculate log likelihood
    > loglikelihood(tree, model)

    -7.5814075725577

### Development
PhyloModels.jl intends to be a full suite of tools for phylogenetic simulation and inference in Julia. Future versions will include
* Further support for heterogenous site rate models
* Tree operators and MCMC methods
