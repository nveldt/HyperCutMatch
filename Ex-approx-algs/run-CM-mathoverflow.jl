using MAT
using SparseArrays
using MatrixNetworks

## Code
include("../src/hypergraph-functions.jl")
include("../include/SparsecardDSFM/hypergraph-clustering-utils.jl")
include("../src/hyper-flow-embed.jl")
include("../src/find-partition.jl")
include("../src/flow-embed-wrappers.jl")

## Load dataset
M = matread("../data/main-data/mathoverflow/math_summary.mat")
topics = M["topics"]
stats = M["stats"]
SetInds = M["SetInds"]
numtimes = 10

for num = 10:length(SetInds)
    mat = matread("../data/main-data/mathoverflow/math_$(SetInds[num])_lcc.mat")
    H = mat["H"]
    m,n = size(H)
    mu = sum(H)
    order = vec(round.(Int64,sum(H,dims = 2)))
    Edges = incidence2elist(H)
    fun = x->((x .> 0).*1)   # all or nothing cut function
    EdgesW = generic_EdgesW(order,fun)

    # For each dataset, we try three different iteration numbers
    # and run the algorithm numtimes times

    Bounds = zeros(3,numtimes)
    Approxs = zeros(3,numtimes)
    Runtimes = zeros(3,numtimes)

    print("$num \t $n \t $m \t $mu \t $(topics[num])")
    T1 = 10*round(Int64,log2(n))
    T2 = 20*round(Int64,log2(n))
    T3 = 30*round(Int64,log2(n))
    Ts = [T1; T2; T3]
    for jj = 1:numtimes
        for it = 1:3
            T = Ts[it]
            tic = time()
            LBs, Lams, RuntimesLocal, Alphas, Approx, H, S = hypergraphcutmatch(Edges,EdgesW,ones(n),T,2,false);
            solvetime = time()-tic

            # double check expansion computation of best set
            expS1 = minimum(Alphas)
            eS = zeros(n)
            eS[S] .= 1
            expS2, sc2 = aon_expansion(Edges,eS) 
            expSr = round(expS2,digits = 3)
            @assert(expS1 == expS2)

            # Recompute the best bound
            ApproxFac, Lbound = lowerbound(H,Alphas,ones(n))
            ratio = round(ApproxFac,digits = 3)
            lb = round(Lbound,digits=3)
            Bounds[it,jj] = Lbound
            Approxs[it,jj] = ApproxFac
            Runtimes[it,jj] = solvetime
            println("$jj \t $it \t $lb \t $expSr \t $ratio \t $solvetime")
        end
    end

    
    matwrite("Output-CM-mathoverflow/math_$(SetInds[num])_output.mat", Dict("Approxs"=>Approxs,"Runtimes"=>Runtimes,"Bounds"=>Bounds))
end

