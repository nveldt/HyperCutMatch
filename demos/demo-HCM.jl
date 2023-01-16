using MAT
using SparseArrays
using MatrixNetworks

## Code
include("../src/HyperCutMatching.jl");

## Load an Amazon hypergraph
M = matread("../data/larger-hypergraphs/Amazon9_lcc.mat")
names = M["names"]
num  = 4    # Choose from 1 to 9
mat = matread("../data/amazon-9/Am_$(num)_lcc.mat")
H = mat["H"]
mu = sum(H)
order = vec(round.(Int64,sum(H,dims = 2)))
Edges = incidence2elist(H)
println("$num \t $n \t $m \t $mu \t hypergraph: $(names[num])")

## What type of hypergraph pi-expansion do you want to approximate?

# Define a cut function: this assumes a symmetric cut function
fun = x->((x .> 0).*1)                  # all-or-nothing cut function
EdgesW = generic_EdgesW(order,fun)

# Choose a node weight function
nodeweights = ones(n);

## Run the Hypergraph Cut Matching algorithm
verbose = true
returnH = true
T = 20                  # number of iterations
eigtol = 1e-8           # tolerance for underlying eigenvalue solver

LBs, Lams, RuntimesLocal, Alphas, Approx, H, S = HyperCutMatch(Edges,EdgesW,nodeweights,T,eigtol,verbose,returnH);

# Check final lower boundm, two different ways
ApproxFactor, LowerBound = lowerbound(H,Alphas,nodeweights)
AppFac = minimum(Alphas)/maximum(LBs)