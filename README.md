# PhyloModels.jl
[![DOI](https://zenodo.org/badge/75206540.svg)](https://zenodo.org/badge/latestdoi/75206540)
[![Latest Release](https://img.shields.io/github/release/jangevaare/PhyloModels.jl.svg)](https://github.com/jangevaare/PhyloModels.jl/releases/latest)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](https://github.com/jangevaare/PhyloModels.jl/blob/master/LICENSE)
[![Build Status](https://travis-ci.org/jangevaare/PhyloModels.jl.svg?branch=master)](https://travis-ci.org/jangevaare/PhyloModels.jl)
[![codecov.io](http://codecov.io/github/jangevaare/PhyloModels.jl/coverage.svg?branch=master)](http://codecov.io/github/jangevaare/PhyloModels.jl?branch=master)

PhyloModels.jl is a package for performing genetic sequence simulation from specified phylogenetic trees (using [PhyloTrees.jl](https://github.com/jangevaare/PhyloTrees.jl) for phylogenetic tree specification, and [SubstitutionModels.jl](https://github.com/BioJulia/SubstitutionModels.jl) for nucleic acid substitution model specification), and for loglikelihood calculation of phylogenetic trees using Felsenstein's tree pruning algorithm¹.

¹ Felsenstein, J. (1981). Evolutionary trees from DNA sequences: a maximum likelihood approach. *Journal of molecular evolution, 17*(6), 368-376.

## Installation

The current release can be installed from the Julia REPL with:

```julia
pkg> add PhyloModels
```

The development version (master branch) can be installed with:

```julia
pkg> add PhyloModels#master
```
