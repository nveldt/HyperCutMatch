using MAT
using SparseArrays
using MatrixNetworks

## Code
include("../src/card-based-CE.jl")

## Run the smart clique expansion method
datasets = ["Newsgroups", "Mushrooms", "Covertype45", "Covertype67"]
for i = 4
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

    alphas = collect(round.(LinRange(0.000,0.04,11),digits = 3))
    for a = 1:length(alphas)
        alpha = alphas[a]
        EdgesW = Alpha_EdgesW(order,alpha)  
        nodeweights = generalized_degree(Edges,EdgesW,n)

        smartsweep = true
        graphnormalize = true
        S, condS, reducetime, sweeptime, y = Smart_CE_Spectral(Edges,EdgesW,nodeweights,smartsweep,graphnormalize)

        matwrite("Output/CE_$(dataset)_alpha_$(alpha)_gnorm_$graphnormalize.mat",Dict("S"=>S,"condS"=>condS,"reducetime"=>reducetime,"sweeptime"=>sweeptime,"y"=>y,"nodeweights"=>nodeweights))
        lS = length(S)
        if lS > n/2
            lS = n-length(S)
        end
        println("$alpha \t $lS $condS $(sweeptime+reducetime)")
    end
end