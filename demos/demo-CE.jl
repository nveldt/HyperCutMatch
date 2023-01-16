using MAT
using SparseArrays
using MatrixNetworks

## Code
include("../src/card-based-CE.jl")

## Load an Amazon hypergraph
num  = 9  
mat = matread("../data/amazon-9/Am_$(num)_lcc.mat")
H = mat["H"]

## Load a different hypergraph
mat = matread("../data/benchmark-hypergraphs/Newsgroups_H.mat")
H = mat["H"]

## Hypergraph statistics
m,n = size(H)
mu = sum(H)
order = vec(round.(Int64,sum(H,dims = 2)))
Edges = incidence2elist(H)
println("$num \t $n \t $m \t $mu \t hypergraph: $(names[num])")

## What type of hypergraph pi-expansion do you want to approximate?
# Define a cut function: this assumes a symmetric cut function

# all-or-nothing cut function
fun = x->((x .> 0).*1)   
EdgesW = generic_EdgesW(order,fun)

# alpha-paramterized function of Li & Milenkovic 
alpha = 0.036
EdgesW = limil_EdgesW(order,alpha)

# Choose a node weight function (e.g., some notion of degrees)
d = vec(sum(H,dims = 1))
dgen = generalized_degree(Edges,EdgesW,n)
nodeweights = d

## Run the smart clique expansion method, normalizing eigenvector by reduced graph degrees
smartsweep = true
graphnormalize = true
S_1, condS_1, reducetime_1, sweeptime_1, y_1 = Smart_CE_Spectral(Edges,EdgesW,nodeweights,smartsweep,graphnormalize)

## Run the smart clique expansion method, normalizing eigenvector by hypergraph degrees
smartsweep = true
graphnormalize = false
S_2, condS_2, reducetime_2, sweeptime_2, y_2 = Smart_CE_Spectral(Edges,EdgesW,nodeweights,smartsweep,graphnormalize)

## Run without the smart sweep feature (often faster, but never better in terms of ratio cut score)
smartsweep = false
S, condS, reducetime, sweeptime, y = Smart_CE_Spectral(Edges,EdgesW,nodeweights,smartsweep)

## Compare against hypergraph cut matching
verbose = true
returnH = true
T = 20                  # number of iterations
eigtol = 1e-4           # tolerance for underlying eigenvalue solver

LBs, Lams, RuntimesLocal, Alphas, Approx, H_exp, S = HyperCutMatch(Edges,EdgesW,nodeweights,T,eigtol,verbose,returnH);

hcmtime = RuntimesLocal[end,3]
cond_hcm = minimum(Alphas)
ApproxFactor, LowerBound = lowerbound(H,Alphas,nodeweights)
