using MAT
using SparseArrays
using MatrixNetworks

## Code
include("../src/HyperCutMatching.jl")

## Run the smart clique expansion method
datasets = ["Newsgroups", "Mushrooms", "Covertype45", "Covertype67"]
for jj = 3
for i = 2
    dataset = datasets[i]
    mat = matread("../data/benchmark-hypergraphs/$(dataset)_H.mat")
    H = mat["H"]
    H, d, order = largest_cc_star(H)
    order = round.(Int64,order)
    mu = sum(order)
    m,n = size(H)
    maxr = round(Int64,maximum(order))
    meanr = round(sum(order)/length(order),digits = 1)
    println("$(dataset[1:4]) \t  $n \t  $m \t  $maxr \t $meanr \t $mu")
    volA = sum(d)
    Edges = incidence2elist(H)

    verbose = false
    returnH = true
    eigtol = 1e-4           # tolerance for underlying eigenvalue solver

    alphas = collect(round.(LinRange(0.000,0.04,11),digits = 3))
    # for a = 1:length(alphas)
    for a = 3
        alpha = alphas[a]

        EdgesW = Alpha_EdgesW(order,alpha)     # hyperedge splitting functions
        nodeweights = generalized_degree(Edges,EdgesW,n)
        T = 10*round(Int64,log2(n))            # number of cut matching iterations
        tic = time()
        LBs, Lams, RuntimesHCM, Conds, Approx, H, S =  HyperCutMatch(Edges,EdgesW,nodeweights,T,eigtol,verbose,returnH);
        solvetime = time()-tic
        ApproxFac, Lbound = lowerbound(H,Conds,nodeweights)

        ratio = round(ApproxFac,digits = 3)
        lb = round(Lbound,digits=3)
        bestcond = round(minimum(Conds),digits = 3)
        numS = min(length(S), n-length(S))
        matwrite("Output/HCM_$(dataset)_alpha_$(alpha)_$(jj).mat",Dict("ApproxFac"=>ApproxFac,"Lbound"=>Lbound,"LBs"=>LBs, "Lams"=>Lams, "RuntimesHCM"=>RuntimesHCM,"Conds"=>Conds,"Approx"=>Approx,"S"=>S,"solvetime"=>solvetime))
        println("\t$alpha: $lb $bestcond $ratio $solvetime $numS $jj")
    end
end
end