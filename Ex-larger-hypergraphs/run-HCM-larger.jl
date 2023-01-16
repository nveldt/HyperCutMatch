using MAT
using SparseArrays
using MatrixNetworks
using LinearAlgebra

## Code
include("../src/HyperCutMatching.jl");

## Load dataset
datasets = ["Cities_H", "Zoo_H", "mathoverflow-answers-all","tripadvisor","Amazon9_lcc","trivago-cities-all"]
numtimes = 5
for jj = 1:numtimes
    for i = 6
        hyper = datasets[i]
        mat = matread("../data/larger-hypergraphs/$hyper.mat")
        H = mat["H"]
        H, d, order = largest_cc_star(H)
        order = round.(Int64,order)
        mu = sum(order)
        m,n = size(H)
        maxr = round(Int64,maximum(order))
        meanr = round(sum(order)/length(order),digits = 1)
        # println("$(hyper[1:4]) \t  $n \t  $m \t  $maxr \t $meanr \t $mu")

        Edges = incidence2elist(H)
        fun = x->((x .> 0).*1)   # all or nothing cut function
        EdgesW = generic_EdgesW(order,fun)

        T = 5*round(Int64,log2(n))

        verbose = true
        returnH = true
        eigtol = 1e-4           # tolerance for underlying eigenvalue solver
        nodeweights = ones(n)
        tic = time()
        LBs, Lams, RuntimesHCM, Alphas, Approx, H, S =  HyperCutMatch(Edges,EdgesW,nodeweights,T,eigtol,verbose,returnH);
        solvetime = time()-tic

        ## double check expansion computation of best set
        expS1 = minimum(Alphas)
        eS = zeros(n)
        eS[S] .= 1
        expS2, sc2 = aon_expansion(Edges,eS)
        @show expS1, expS2
        @assert(expS1 == expS2)

        # Recompute the best bound
        println("a sparse eigevec comp")
        ApproxFac, Lbound = lowerbound(H,Alphas,ones(n))
        matwrite("Output/HCM-$(hyper)-$(jj).mat",Dict("ApproxFac"=>ApproxFac,"Lbound"=>Lbound,"LBs"=>LBs, "Lams"=>Lams, "RuntimesHCM"=>RuntimesHCM, "Alphas"=>Alphas,"Approx"=>Approx,"S"=>S))

        ratio = round(ApproxFac,digits = 3)
        lb = round(Lbound,digits=3)
        println("$(hyper[1:4]) \t $lb \t $ratio \t $solvetime $jj")

    end
end



##
hyper = "trivago-cities-all"


# mat = matread("../data/larger-hypergraphs/$hyper.mat")
# H = mat["H"]
# H, d, order = largest_cc_star(H)
# order = round.(Int64,order)
# mu = sum(order)
# m,n = size(H)
# maxr = round(Int64,maximum(order))
# meanr = round(sum(order)/length(order),digits = 1)
# # println("$(hyper[1:4]) \t  $n \t  $m \t  $maxr \t $meanr \t $mu")

# Edges = incidence2elist(H)
# fun = x->((x .> 0).*1)   # all or nothing cut function
# EdgesW = generic_EdgesW(order,fun)

# include("../src/card-based-CE.jl")
# tic = time()
# A,p = MinDistortionCliqueExpansion(Edges,EdgesW,n)
# @show time()-tic