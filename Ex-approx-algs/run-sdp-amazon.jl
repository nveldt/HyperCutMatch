using MAT
using SparseArrays
using MatrixNetworks

## Code
include("../src/hypergraph-helper-functions.jl")
include("../src/hypergraph-sc-sdp.jl")
include("../src/hypergraph-sc-lp.jl")

## Load datasets and run

# WARNING: Calling Mosek from Julia/JuMP is worse than using CVX in Matlab.
# CVX sets up and solves the dual problem. This is very costly, but 
# without this step, the problem crashes. 

M = matread("../data/larger-hypergraphs/Amazon9_lcc.mat")
names = M["names"]

for num = 3:5
    mat = matread("../data/amazon-9/Am_$(num)_lcc.mat")
    H = mat["H"]
    m,n = size(H)
    order = round.(Int64,vec(sum(H,dims=2)))
    d = round.(Int64,vec(sum(H,dims=1)))
    Elist = incidence2elist(H)

    print("$num \t $n \t $m ")

    ## Run SDP solver
    outputflag = true
    time_limit = 1800.0
    optimal, X, Y, solverstats = HypergraphExpansionSDP(H,time_limit,outputflag)
    if optimal
        obj = solverstats[1]
        lb = round(obj/2,digits = 3)
        setuptime = round(solverstats[2],digits = 3)
        solvetime = round(solverstats[3],digits = 3)
        c, expS = round_hyperLR(X,Elist)
        expS = round(expS,digits = 3)
        ratio = round(expS/lb,digits = 3)
        println(" \t $lb \t $expS \t $ratio \t $solvetime")
    else
        setuptime = round(solverstats[1],digits = 3)
        solvetime = round(solverstats[2],digits = 3)
        c = 0
        expS = 0
        lb = 0
        println(" \t timed out is: $X \t $solvetime")
    end
    matwrite("Output/amazon/Am_$(num)_SDP_solution_primal.mat", Dict("X"=>X,"Y"=>Y, "c"=>c, "expS"=>expS,"lb"=>lb, "optimal"=>optimal,"solverstats"=>solverstats))
end

