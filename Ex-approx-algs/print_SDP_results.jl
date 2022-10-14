using MAT

include("../src/hypergraph-sc-lp.jl")


## Output results for mathoverflow
M = matread("../data/main-data/mathoverflow/math_summary.mat")
topics = M["topics"]
stats = M["stats"]
SetInds = M["SetInds"]

# Need to report manually the solution/times for SDP solver
times = zeros(length(SetInds))
times[1] = 20.10
times[5] = 82.28  
solved = zeros(Bool,length(SetInds))
solved[1] = true
solved[5] = true

for num = 1:5
    tpc = topics[num]
    mat = matread("../data/main-data/mathoverflow/math_$(SetInds[num])_lcc.mat")
    
    H = mat["H"]
    Elist = incidence2elist(H)
    m,n = size(H)
    order = round.(Int64,vec(sum(H,dims=2)))
    mu = sum(order) 
    print("$num \t $n \t $m \t $mu \t ")
            
    mat2 = matread("../data/main-data/mathoverflow/math_$(SetInds[num])_SDP_solution.mat")
    X = mat2["X"]
    X = X+ X'
    totaltime = mat2["timeSDP"]
    setuptime = totaltime - times[num]
    lb = mat2["lowerbound"]
    if solved[num]
        eS, expS = round_hyperLR(X,Elist)
        expSr = round(expS,digits = 3)
        lbr = round(lb,digits = 3)
        approx = round(expS/lb,digits = 3)
        println("| $lbr \t $expSr \t $approx \t $(times[num])")
    else
        println("Failed")
    end
    
end

## Output results for trivago
M = matread("../data/main-data/trivago/trivago_summary.mat")
topics = M["cities"]
stats = M["stats"]
SetInds = M["SetInds"]

# Manually report the solution/times for SDP solver from MATLAB OUTPUT
times = zeros(length(SetInds))
solved = zeros(Bool,length(SetInds))
times[1] = 23.15  
solved[1] = true

times[8] = 212.39 
solved[8] = true  

times[10] = 280.31  # set 2718
solved[10] = true  

for num = 1:12
    tpc = topics[num]
    mat = matread("../data/main-data/trivago/trivago_$(SetInds[num])_lcc.mat")
    
    H = mat["H"]
    Elist = incidence2elist(H)
    m,n = size(H)
    order = round.(Int64,vec(sum(H,dims=2)))
    mu = sum(order) 
    print("$num \t $n \t $m \t $mu \t ")
            
    mat2 = matread("../data/main-data/trivago/trivago_$(SetInds[num])_SDP_solution.mat")
    X = mat2["X"]
    X = X+ X'
    totaltime = mat2["timeSDP"]
    setuptime = totaltime - times[num]
    lb = mat2["lowerbound"]
    if solved[num]
        eS, expS = round_hyperLR(X,Elist)
        expSr = round(expS,digits = 3)
        lbr = round(lb,digits = 3)

        approx = round(expS/lb,digits = 3)
        println("| $lbr \t $expSr \t $approx \t $(times[num]) \t $setuptime")
    else
        println("Failed")
    end
    
end