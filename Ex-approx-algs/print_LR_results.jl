using MAT

include("../src/hypergraph-sc-lp.jl")

M = matread("../data/main-data/amazon/Amazon9_lcc.mat")
names = M["names"]

println("Amazon hypergraph results")
for num = 1:5
mat = matread("../data/main-data/amazon/Am_$(num)_lcc.mat")
H = mat["H"]
Elist = incidence2elist(H)
m,n = size(H)
order = round.(Int64,vec(sum(H,dims=2)))
mu = sum(order) 
print("$num \t $n \t $m \t $mu \t ")
        
mat2 = matread("Output-LR-Amazon/amazon_$(num)_output.mat")
X = mat2["X"]
X = X+X'
lb = mat2["solverstats"][1]/2
optimal = mat2["optimal"]
solverstats = mat2["solverstats"]
if optimal
    m = 10
    eS, expS = round_LP_many(X,m)
    expSr = round(expS,digits = 3)
    # @show expS, lb
    approx = round(expS/lb,digits = 3)
    setuptime = round(solverstats[2],digits = 3)
    solvetime = round(solverstats[3],digits = 3)
    println("| $lb \t $expSr \t $approx \t $solvetime")
else
    println("Failed")
end

end

## Output results for mathoverflow
M = matread("../data/main-data/trivago/trivago_summary.mat")
topics = M["cities"]
stats = M["stats"]
SetInds = M["SetInds"]

for num = 1:length(SetInds)
    tpc = topics[num]
    mat = matread("../data/main-data/trivago/trivago_$(SetInds[num])_lcc.mat")
    
    H = mat["H"]
    Elist = incidence2elist(H)
    m,n = size(H)
    order = round.(Int64,vec(sum(H,dims=2)))
    mu = sum(order) 
    print("$num \t $n \t $m \t $mu \t ")
            
    mat2 = matread("Output-LR-Trivago/trivago_$(SetInds[num])_output.mat")
    X = mat2["X"]
    X = X+X'
    lb = mat2["lb"]
    optimal = mat2["optimal"]
    solverstats = mat2["solverstats"]
    ebad = mat2["expS"]
    approxbad = round(ebad/lb,digits = 3)
    if optimal
        m = 10
        eS, expS = round_LP_many(X,m)
        expSr = round(expS,digits = 3)
        approx = round(expS/lb,digits = 3)
        approxrep = min(approx,approxbad) # report best of both
        setuptime = round(solverstats[2],digits = 3)
        solvetime = round(solverstats[3],digits = 3)
        println("| $lb \t $expSr \t $approxrep \t $approxbad \t $solvetime \t $tpc")
    else
        println("Failed \t \t \t $tpc")
    end
    
end


## Output results for mathoverflow
M = matread("../data/main-data/mathoverflow/math_summary.mat")
topics = M["topics"]
stats = M["stats"]
SetInds = M["SetInds"]

for num = 1:length(SetInds)
    tpc = topics[num]
    mat = matread("../data/main-data/mathoverflow/math_$(SetInds[num])_lcc.mat")
    
    H = mat["H"]
    Elist = incidence2elist(H)
    m,n = size(H)
    order = round.(Int64,vec(sum(H,dims=2)))
    mu = sum(order) 
    print("$num \t $n \t $m \t $mu \t ")
            
    mat2 = matread("Output-LR-mathoverflow/math_$(SetInds[num])_output.mat")
    X = mat2["X"]
    X = X+X'
    lb = mat2["lb"]
    optimal = mat2["optimal"]
    solverstats = mat2["solverstats"]
    ebad = mat2["expS"]
    approxbad = round(ebad/lb,digits = 3)
    if optimal
        m = 10
        eS, expS = round_LP_many(X,m)
        expSr = round(expS,digits = 3)
        approx = round(expS/lb,digits = 3)
        setuptime = round(solverstats[2],digits = 3)
        solvetime = round(solverstats[3],digits = 3)
        println("| $lb \t $expSr \t $approx \t $approxbad \t $solvetime \t $tpc")
    else
        println("Failed \t \t \t $tpc")
    end
    
end