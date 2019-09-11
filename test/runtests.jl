using Test,
      PhyloModels

# An example from "Molecular Evolution: A Statistical Approach" by Ziheng Yang

# Describe a Phylogenetic Tree
tree = Tree()

# Add 9 nodes
addnodes!(tree, 9)

# Connect with branches...
addbranch!(tree, 9, 6, 0.1)
addbranch!(tree, 9, 8, 0.1)
addbranch!(tree, 6, 7, 0.1)
addbranch!(tree, 6, 3, 0.2)
addbranch!(tree, 7, 1, 0.2)
addbranch!(tree, 7, 2, 0.2)
addbranch!(tree, 8, 4, 0.2)
addbranch!(tree, 8, 5, 0.2)

# Specify sequences of leaf nodes
node_data = NodeDNA()
node_data[1] = DNASeq("T")
node_data[2] = DNASeq("C")
node_data[3] = DNASeq("A")
node_data[4] = DNASeq("C")
node_data[5] = DNASeq("C")

# Specify a Nucleic Acid Substitution Model
model = K80(2.0)

# loglikelihood calculation
ll = loglikelihood(tree, model, node_data)

@test ll == -7.5814075725577
