using MAT
using SparseArrays
using MatrixNetworks
using LinearAlgebra

## Code
include("../src/hypergraph-functions.jl")
include("../include/SparsecardDSFM/hypergraph-clustering-utils.jl")
include("../src/hyper-flow-embed.jl")
include("../src/find-partition.jl")
include("../src/flow-embed-wrappers.jl")
include("../src/hypergraph-sc-lp.jl")

## Load dataset
datasets = ["Cities_H", "Zoo_H", "mathoverflow-answers-all","tripadvisor","Amazon9_lcc","trivago-cities-all"]

numtimes = 5
for jj = 1:numtimes
    for i = 1:5
        hyper = datasets[i]
        mat = matread("../data/larger/$hyper.mat")
        H = mat["H"]
        H, d, order = largest_cc_star(H)
        order = round.(Int64,order)
        mu = sum(order)
        m,n = size(H)
        maxr = round(Int64,maximum(order))
        meanr = round(sum(order)/length(order),digits = 1)
        println("$(hyper[1:4]) \t  $n \t  $m \t  $maxr \t $meanr \t $mu")

        Edges = incidence2elist(H)
        fun = x->((x .> 0).*1)   # all or nothing cut function
        EdgesW = generic_EdgesW(order,fun)

        T1 = 5*round(Int64,log2(n))
        T3 = 25*round(Int64,log2(n))
        Ts = [T1; T3]

        for it = 2
            eigtol = 1e-4
            T = Ts[it]
            @show T, it
            nodeweights = ones(n)
            tic = time()
            LBs, Lams, RuntimesHCM, Alphas, Approx, H, S = hypergraphcutmatch(Edges,EdgesW,nodeweights,T,2,eigtol,true);
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
            ApproxFac, Lbound = lowerbound(H,Alphas,ones(n),false)
            matwrite("Output/HCM-$(hyper)-$(it)-$(jj).mat",Dict("ApproxFac"=>ApproxFac,"Lbound"=>Lbound,"LBs"=>LBs, "Lams"=>Lams, "RuntimesHCM"=>RuntimesHCM, "Alphas"=>Alphas,"Approx"=>Approx,"H"=>H,"S"=>S))

            ratio = round(ApproxFac,digits = 3)
            lb = round(Lbound,digits=3)
            # println("$lb \t $expSr \t $ratio \t $solvetime")
        end
    end
end