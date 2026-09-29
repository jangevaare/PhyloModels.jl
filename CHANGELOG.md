# Changelog

This history was reconstructed from the commits and source changes between
successive tags. Changes after the latest tag are listed as unreleased, even
where the package version has already been updated.

## Unreleased (version set to 0.4.0)

- Reuse likelihood calculations through a `LikelihoodWorkspace` when evaluating
  the same tree repeatedly. Identical alignment columns are grouped by default,
  while returned per-site calculations still follow the original alignment.
- Reuse transition probabilities during sequence simulation, especially for
  repeated site-rate categories and common nucleotide states.
- Accept an explicit random number generator for sequence generation and
  simulation, so callers can control the random stream.
- Allow newer GeneticBitArrays, PhyloTrees, SubstitutionModels, and StatsBase
  releases. Refresh the test and dependency-update workflows.

## [0.3.4](https://github.com/jangevaare/PhyloModels.jl/compare/v0.3.3...v0.3.4)

- Allow PhyloTrees 0.11 and SubstitutionModels 0.5 alongside earlier supported
  versions, and correct the Julia compatibility declaration.
- Move package tests from Travis CI to GitHub Actions and run them against
  stable, long-term-support, and nightly Julia releases.

## [0.3.3](https://github.com/jangevaare/PhyloModels.jl/compare/v0.3.2...v0.3.3)

- Revert a restrictive node-data type annotation added in 0.3.2, restoring
  likelihood calculations for the supported sequence dictionaries.
- Add a DOI badge to the README.

## [0.3.2](https://github.com/jangevaare/PhyloModels.jl/compare/v0.3.1...v0.3.2)

- Update the likelihood function's node-data type annotation. This change was
  reversed in 0.3.3.
- Add code coverage reporting, Julia 1.4 testing, and clearer installation
  instructions.

## [0.3.1](https://github.com/jangevaare/PhyloModels.jl/compare/v0.3.0...v0.3.1)

- Fix the exported spelling of `simulate!`, so the function is available
  through `using PhyloModels`.
- Update compatibility bounds for SubstitutionModels, PhyloTrees, and
  GeneticBitArrays. Add automated dependency and release workflows.

## [0.3.0](https://github.com/jangevaare/PhyloModels.jl/compare/v0.2.2...v0.3.0)

- Rebuild for Julia 1.x using `Project.toml`.
- Move nucleotide substitution models into SubstitutionModels.jl, tree structures
  into PhyloTrees.jl, and sequence storage into GeneticBitArrays.jl. Re-export
  their public names from PhyloModels.
- Adapt simulation and tree likelihood calculations to the new sequence and
  model types. Leaf sequences are supplied separately from the tree; likelihood
  calculations can optionally be returned with the result.

## [0.2.2](https://github.com/jangevaare/PhyloModels.jl/compare/v0.2.1...v0.2.2)

- Correct the GTR rate matrix and the order of matrix multiplication used when
  simulating descendant sequences.
- Introduce a GTR transition-probability method for an array of times. As
  tagged, this method still refers to an undefined variable and cannot be
  used successfully.
- Remove the experimental sandbox API and update CI configuration.

## [0.2.1](https://github.com/jangevaare/PhyloModels.jl/compare/v0.2.0...v0.2.1)

- Reduce work in tree likelihood calculations.
- Add indexing support for experimental trace objects and `deleteat!` for
  their stored iterations.

## [0.2.0](https://github.com/jangevaare/PhyloModels.jl/compare/v0.1.0...v0.2.0)

- Adapt simulation and likelihood calculations to changes in PhyloTrees:
  observed sequences now live in a node-indexed dictionary instead of inside
  tree nodes.
- Add `simulate(n, model)` to draw a root sequence from model frequencies,
  and add copying methods for sequences.
- Check for missing leaf observations and multiple roots before calculating a
  tree likelihood. Fix sequence simulation and update the examples.

## [0.1.0](https://github.com/jangevaare/PhyloModels.jl/tree/v0.1.0)

- Initial release with DNA sequence simulation and tree likelihood calculation.
- Include JC69, K80, F81, F84, HKY85, TN93, and GTR substitution models, plus
  experimental prior and proposal code.
