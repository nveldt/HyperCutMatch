using MAT
using SparseArrays
using MatrixNetworks
using LinearAlgebra

## Code
include("../src/HyperCutMatching.jl");

## Load dataset
M = matread("../data/mathoverflow-small/math_summary.mat")
topics = M["topics"]
stats = M["stats"]
SetInds = M["SetInds"]
numtimes = 10

for num = 1:length(SetInds)
    mat = matread("../data/mathoverflow-small/math_$(SetInds[num])_lcc.mat")
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
    T2 = 30*round(Int64,log2(n))
    Ts = [T1; T2]

    verbose = false
    returnH = true
    eigtol = 1e-8           # tolerance for underlying eigenvalue solver

    for jj = 1:numtimes
        for it = 1:2
            T = Ts[it]
            tic = time()
            LBs, Lams, RuntimesLocal, Alphas, Approx, H, S = HyperCutMatch(Edges,EdgesW,ones(n),T,eigtol,verbose,returnH);
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

    
    matwrite("Output/mathoverflow/math_$(SetInds[num])_HCM_solution.mat", Dict("Approxs"=>Approxs,"Runtimes"=>Runtimes,"Bounds"=>Bounds))
end

