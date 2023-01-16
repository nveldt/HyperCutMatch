using MAT
using SparseArrays
using MatrixNetworks
using LinearAlgebra

## Code
include("../src/HyperCutMatching.jl");
numtimes = 10

## Load an Amazon hypergraph
M = matread("../data/larger-hypergraphs/Amazon9_lcc.mat")
names = M["names"]

for num = 1:9

    mat = matread("../data/amazon-9/Am_$(num)_lcc.mat")
    H = mat["H"]
    mu = sum(H)
    m,n = size(H)
    order = vec(round.(Int64,sum(H,dims = 2)))
    Edges = incidence2elist(H)
    println("$num \t $n \t $m \t $mu \t hypergraph: $(names[num])")
    fun = x->((x .> 0).*1)                 
    EdgesW = generic_EdgesW(order,fun)
    nodeweights = ones(n);

    Bounds = zeros(2,numtimes)
    Approxs = zeros(2,numtimes)
    Runtimes = zeros(2,numtimes)
    T1 = 10*round(Int64,log2(n))
    T2 = 30*round(Int64,log2(n))
    Ts = [T1; T2]

    verbose = false
    returnH = true
    eigtol = 1e-8           # tolerance for underlying eigenvalue solver

    for jj = 1:numtimes
        for it = 1:length(Ts)
            T = Ts[it]
            tic = time()
            LBs, Lams, RuntimesLocal, Alphas, Approx, H, S = HyperCutMatch(Edges,EdgesW,nodeweights,T,eigtol,verbose,returnH);
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

    matwrite("Output/amazon/Am_$(num)_HCM_solution.mat", Dict("Approxs"=>Approxs,"Runtimes"=>Runtimes,"Bounds"=>Bounds))

end