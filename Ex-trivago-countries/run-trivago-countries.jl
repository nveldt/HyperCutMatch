using MAT
using SparseArrays
using MatrixNetworks
using LinearAlgebra

# Code
include("../src/hypergraph-functions.jl")
include("../include/SparsecardDSFM/hypergraph-clustering-utils.jl")
include("../src/hyper-flow-embed.jl")
include("../src/find-partition.jl")
include("../src/flow-embed-wrappers.jl")
include("../src/hypergraph-sc-lp.jl")

# Load smaller dataset
M = matread("../data/trivago-large/trivago_countries_large_summary.mat")
Labels = M["SpecialLabels"]
LabelNames = M["LabelNames"]

inds = [5, 8, 17, 20]
for jj = 5
for i in inds[1:2]
countryind = Labels[i]
countryname = LabelNames[countryind]
l = countryind
mat = matread("../data/trivago-large/trivago_countries_$(l)_lcc.mat")
H = mat["H"] 
k = 2
H, d, order = hyperkcore_multi(H,k)
matwrite("../data/trivago-large/trivago_countries_$(l)_2core.mat", Dict("H"=>H))

order = round.(Int64,order)
mu = sum(order)
m,n = size(H)
maxr = round(Int64,maximum(order))
meanr = round(sum(order)/length(order),digits = 1)
println("$countryname: \t  $n \t  $m \t  $maxr \t $meanr \t $mu")
volA = sum(d)
Edges = incidence2elist(H);

# Run Hypergraph Cut Matching
deltas = [1 1.1 1.5 2 5 10 100]
for a = 1:length(deltas)
    delta = deltas[a]
    EdgesW = Delta_EdgesW(order,delta)     # hyperedge splitting functions
    nodeweights = d                        # node weights (degrees, for conductance)
    T = 10*round(Int64,log2(n))             # number of cut matching iterations
    eigtol = 1e-4                          # tolerance for eigensolver
    eigsolver = 2                          # eigsolver is KrylovKit
    verbose = false                         # silence output 
    tic = time()
    LBs, Lams, RuntimesHCM, Conds, Approx, H, S = hypergraphcutmatch(Edges,EdgesW,nodeweights,T,2,eigtol,verbose);
    solvetime = time()-tic
    ApproxFac, Lbound = lowerbound(H,Conds,nodeweights,false)

    ratio = round(ApproxFac,digits = 3)
    lb = round(Lbound,digits=3)
    bestcond = round(minimum(Conds),digits = 3)
    numS = min(length(S), n-length(S))
    matwrite("Output/HCM_tric_$(l)_2core_$(delta)_$(jj)_longer.mat",Dict("ApproxFac"=>ApproxFac,"Lbound"=>Lbound,"LBs"=>LBs, "Lams"=>Lams, "RuntimesHCM"=>RuntimesHCM,"Conds"=>Conds,"Approx"=>Approx,"H"=>H,"S"=>S,"solvetime"=>solvetime))
    println("\t$delta: $lb $bestcond $ratio $solvetime $numS")
end

end
end