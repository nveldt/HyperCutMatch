using MAT
using SparseArrays
using MatrixNetworks

## Code
include("../src/hypergraph-helper-functions.jl")
include("../src/hypergraph-sc-lp.jl")

## Load datasets and run
M = matread("../data/larger-hypergraphs/Amazon9_lcc.mat")
names = M["names"]

for num = 1:5
    mat = matread("../data/amazon-9/Am_$(num)_lcc.mat")
    H = mat["H"]
    m,n = size(H)
    order = round.(Int64,vec(sum(H,dims=2)))
    d = round.(Int64,vec(sum(H,dims=1)))
    Elist = incidence2elist(H)

    print("$num \t $n \t $m ")
    # Run LP solver
    outputflag = false
    time_limit = 1800.0
    optimal, X, Y, solverstats = HypergraphLeightonRao(H,time_limit,outputflag)
    X=  triu(X)
    X = X+X'
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
    matwrite("Output/amazon/Am_$(num)_LPjump_solution.mat", Dict("X"=>X,"Y"=>Y, "c"=>c, "expS"=>expS,"lb"=>lb, "optimal"=>optimal,"solverstats"=>solverstats))
end


for num = 1:5
    mat = matread("../data/amazon-9/Am_$(num)_lcc.mat")
    H = mat["H"]
    m,n = size(H)
    order = round.(Int64,vec(sum(H,dims=2)))
    d = round.(Int64,vec(sum(H,dims=1)))
    Elist = incidence2elist(H)

    print("$num \t $n \t $m ")
    # Run LP solver
    outputflag = false
    time_limit = 1800.0
    optimal, X, Y, solverstats = HypergraphLR_noJuMP(H,time_limit, outputflag)
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
    matwrite("Output/amazon/Am_$(num)_LP_solution.mat", Dict("X"=>X,"Y"=>Y, "c"=>c, "expS"=>expS,"lb"=>lb, "optimal"=>optimal,"solverstats"=>solverstats))
end
