using MAT
using SparseArrays
using MatrixNetworks

## Code
include("../src/card-based-CE.jl")

## Load a different hypergraph
mat = matread("../data/benchmark-hypergraphs/Newsgroups_H.mat")
H = mat["H"]

## Hypergraph statistics
m,n = size(H)
mu = sum(H)
order = vec(round.(Int64,sum(H,dims = 2)))
Edges = incidence2elist(H)
println("$num \t $n \t $m \t $mu \t hypergraph: $(names[num])")

# alpha-paramterized function of Li & Milenkovic 
alpha = 0.012
EdgesW = limil_EdgesW(order,alpha)

# Choose a node weight function (e.g., some notion of degrees)
d = vec(sum(H,dims = 1))
nodeweights = d

## Run the smart clique expansion method 
smartsweep = true
S, condS_sm, reducetime_sm, sweeptime_sm = Smart_CE_Spectral(Edges,EdgesW,nodeweights,smartsweep)

## Check the sweep cut procedure for a given y
M = matread("checky.mat")
y2 = vec(M["y"])

## Check the eigenvalue computation
A,p = MinDistortionCliqueExpansion(Edges, EdgesW,n)
dA = vec(sum(A,dims = 2))
Dhalf = Diagonal(dA.^(-1/2))
L = I - Dhalf*A*Dhalf
An = Dhalf*A*Dhalf

## Eigenvalue and vec time
sc = n
tic = time()
Vl,Vc,convinfo = eigsolve(L + sc*LinearAlgebra.I, 2, :SR; tol = 1e-8, maxiter = 1000, verbosity = 0)
Ltime = time()-tic

lam2 = Real(Vl[2])-sc
v = Real.(Vc[2])
@show sum(L*v - lam2*v)

## For A
sc = n
tic = time()
Vl,Vc,convinfo = eigsolve(A, 2, :LR; tol = 1e-8, maxiter = 1000, verbosity = 0)
Atime = time()-tic

lam2a = Real(Vl[2])
va = Real.(Vc[2])
@show sum(A*va - lam2a*va)
