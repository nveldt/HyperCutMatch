using MAT
using SparseArrays
using MatrixNetworks

## Code
include("../src/hypergraph-helper-functions.jl")
include("../src/hypergraph-sc-lp.jl")

## Load dataset
M = matread("../data/trivago-small/trivago_summary.mat")
cities = M["cities"]
stats = M["stats"]
SetInds = M["SetInds"]

for num = 36:length(SetInds)


mat = matread("../data/trivago-small/trivago_$(SetInds[num])_lcc.mat")

    H = mat["H"]
    m,n = size(H)
    order = round.(Int64,vec(sum(H,dims=2)))
    d = round.(Int64,vec(sum(H,dims=1)))
    Elist = incidence2elist(H)

    nm = cities[num][1:9]
    print("$num \t $nm \t $n \t $m ")
    # Run LP solver
    outputflag = false
    time_limit = 180000.0
    optimal, X, Y, solverstats = HypergraphLR_noJuMP(H,time_limit,outputflag)
    if optimal
        obj = solverstats[1]
        lb = round(obj/2,digits = 3)
        setuptime = round(solverstats[2],digits = 3)
        solvetime = round(solverstats[3],digits = 3)

        m = 50
        c, expS = round_LP_many(X,m,Elist)
        expS = round(expS,digits = 3)
        ratio = round(expS/lb,digits = 3)
        println(" \t $lb \t $expS \t $ratio \t $solvetime")
    else
        setuptime = round(solverstats[1],digits = 3)
        solvetime = round(solverstats[2],digits = 3)
        c = 0
        expS = 0
        lb = 0
        println(" \t $X \t $solvetime")
        X = 0
    end
    matwrite("Output/trivago/trivago_$(SetInds[num])_LP_solution_notimelimit.mat", Dict("X"=>X,"Y"=>Y, "c"=>c, "expS"=>expS,"lb"=>lb, "optimal"=>optimal,"solverstats"=>solverstats))
end