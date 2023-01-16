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
alpha = 0.012
EdgesW = limil_EdgesW(order,alpha)

# Choose a node weight function (conductance here)
d = vec(sum(H,dims = 1))
dgen = generalized_degree(Edges,EdgesW,n)
nodeweights = dgen
@show norm(d - dgen)
nodeweights = d

## Run the smart clique expansion method 
smartsweep = true
S, condS_sm, reducetime_sm, sweeptime_sm, A, y = Smart_CE_Spectral(Edges,EdgesW,nodeweights,smartsweep)

@show length(S), condS_sm, sweeptime_sm


## Check the sweep cut procedure for a given y
M = matread("checky.mat")
y2 = vec(M["y"])

volA = sum(nodeweights)
S2, ratioS2 = smart_sweep(Edges,EdgesW,nodeweights,y,volA)
S3, ratioS3 = smart_sweep(Edges,EdgesW,nodeweights,y2,volA)

## Check the eigenvalue computation
L,Dhalf = nLaplacian(A,true)
n = size(L,1)
sc = n
@time Vl,Vc,convinfo = eigsolve(L + sc*LinearAlgebra.I, 2, :SR; tol = 1e-8, maxiter = 1000, verbosity = 0)
lam2 = Real(Vl[2])-sc
v = Real.(Vc[2])
@show sum(L*v - lam2*v)

## Run without the smart sweep feature (often faster, but never better in terms of ratio cut score)
smartsweep = false
S, condS, reducetime, sweeptime = Smart_CE_Spectral(Edges,EdgesW,nodeweights,smartsweep)

@show length(S), condS, sweeptime


## Compare against hypergraph cut matching
verbose = true
returnH = true
T = 20                  # number of iterations
eigtol = 1e-4           # tolerance for underlying eigenvalue solver

LBs, Lams, RuntimesLocal, Alphas, Approx, H_exp, S = HyperCutMatch(Edges,EdgesW,nodeweights,T,eigtol,verbose,returnH);

hcmtime = RuntimesLocal[end,3]
cond_hcm = minimum(Alphas)
ApproxFactor, LowerBound = lowerbound(H,Alphas,nodeweights)
