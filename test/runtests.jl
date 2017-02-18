using PhyloTrees
using PhyloModels
using Distributions
using Base.Test

# 1.0 Substitution models
# 1.1 JC69
ta = rand(Gamma(10.))
a = JC69([1.])
@test maximum(abs(expm(Q(a) * ta) .- P(a, ta))) < 1e-14
@test rand(Multinomial(1, a.π))' * P(a, Inf) == a.π'

# 1.2 K80
@test Q(K80([1., 1.])) == Q(a)
tb = rand(Gamma(10.))
b = K80([1., 2.])
@test maximum(abs(expm(Q(b) * tb) .- P(b, tb))) < 1e-14
@test rand(Multinomial(1, b.π))' * P(b, Inf) == b.π'

# 1.3 F81
@test Q(F81([4.], [0.25, 0.25, 0.25, 0.25])) == Q(a)
πc = [0.3, 0.3, 0.2, 0.2]
tc = rand(Gamma(10.))
c = F81([1.], πc)
@test maximum(abs(expm(Q(c) * tc) .- P(c, tc))) < 1e-14
@test rand(Multinomial(1, c.π))' * P(c, Inf) == c.π'

# 1.4 F84
πd = [0.3, 0.3, 0.2, 0.2]
td = rand(Gamma(10.))
d = F84([1., 2.], πd)
@test maximum(abs(expm(Q(d) * td) .- P(d, td))) < 1e-14
@test rand(Multinomial(1, d.π))' * P(d, Inf) == d.π'

# 1.5 HKY85
πe = [0.3, 0.3, 0.2, 0.2]
te = rand(Gamma(10.))
e = HKY85([1., 2.], πe)
@test maximum(abs(expm(Q(e) * te) .- P(e, te))) < 1e-14
@test rand(Multinomial(1, e.π))' * P(e, Inf) == e.π'

# 1.6 TN93
πf = [0.3, 0.3, 0.2, 0.2]
tf = rand(Gamma(10.))
f = TN93([1., 2., 3.], πf)
@test maximum(abs(expm(Q(f) * tf) .- P(f, tf))) < 1e-14
@test rand(Multinomial(1, f.π))' * P(f, Inf) == f.π'

# 2.0 Simulation
g = Tree()
addnode!(g)
branch!(g, 1, 10.0)
branch!(g, 1, 5.0)
branch!(g, 2, 20.0)

model = JC69([1.0e-5])
node_data = Dict{Int64, Sequence}()
node_data[1] = simulate(1000, model)
simulate!(node_data, g, JC69([1.0e-5]), rand(Gamma(1.), 1000))
for i = 1:length(node_data)
  @test length(node_data[i]) == 1000
end

# 3.0 Inference
# 3.1 Tree log likelihood

# Test from Section 4.2 of
# Molecular Evolution: A Statistical Approach, Ziheng Yang

# Build tree
tree = Tree()
addnodes!(tree, 9)
addbranch!(tree, 9, 6, 0.1)
addbranch!(tree, 9, 8, 0.1)
addbranch!(tree, 6, 7, 0.1)
addbranch!(tree, 6, 3, 0.2)
addbranch!(tree, 7, 1, 0.2)
addbranch!(tree, 7, 2, 0.2)
addbranch!(tree, 8, 4, 0.2)
addbranch!(tree, 8, 5, 0.2)

# Set state of leaf nodes
node_data = Dict{Int64, Sequence}()
node_data[1] = Sequence("T")
node_data[2] = Sequence("C")
node_data[3] = Sequence("A")
node_data[4] = Sequence("C")
node_data[5] = Sequence("C")

# Parametrize substitution model
model = K80([2.])

# Calculate log likelihood
@test loglikelihood(tree, model, node_data) == -7.5814075725577
