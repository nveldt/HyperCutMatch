using MAT
using StatsBase
using LaTeXStrings

include("../src/hypergraph-sc-lp.jl")

M = matread("../data/main-data/mathoverflow/math_summary.mat")
topics = M["topics"]
stats = M["stats"]
SetInds = M["SetInds"]
numtimes = 10

# Manually report of how the SDP solver did
times = zeros(length(SetInds))
times[1] = 20.10
times[5] = 82.28  
solved = zeros(Bool,length(SetInds))
solved[1] = true
solved[5] = true 


approx_cm1 = Vector{Float64}()
runs_cm1 = Vector{Float64}()

approx_cm2 = Vector{Float64}()
runs_cm2 = Vector{Float64}()

approx_cm3 = Vector{Float64}()
runs_cm3 = Vector{Float64}()

approx_lp = Vector{Float64}()
runs_lp = Vector{Float64}()

approx_sdp = Vector{Float64}()
runs_sdp = Vector{Float64}()

N = Vector{Int64}()
M = Vector{Int64}()
Mu = Vector{Int64}()

for num = 1:length(SetInds)

    # Load the data
    mat = matread("../data/main-data/mathoverflow/math_$(SetInds[num])_lcc.mat")
    H = mat["H"]
    m,n = size(H)
    mu = round(Int64,sum(H))
    Elist = incidence2elist(H)
    push!(N,n)
    push!(M,m)
    push!(Mu,mu)

    # Load results from cut-matching
    mat_cm = matread("Output-CM-mathoverflow/math_$(SetInds[num])_output.mat")
    Approx = mat_cm["Approxs"]
    Runtimes = mat_cm["Runtimes"]
    Bounds = mat_cm["Bounds"]

    push!(approx_cm1,mean(Approx[1,:]))
    push!(approx_cm2,mean(Approx[2,:]))
    push!(approx_cm3,mean(Approx[3,:]))
    @assert(length(Approx[1,:]) == 10)

    push!(runs_cm1,mean(Runtimes[1,:]))
    push!(runs_cm2,mean(Runtimes[2,:]))
    push!(runs_cm3,mean(Runtimes[3,:]))

    # Load results from LP relaxation
    try 
        mat_cm = matread("Output-LR-mathoverflow/math_$(SetInds[num])_output.mat")
        X = mat_cm["X"]
        X = X+X'
        lb = mat_cm["lb"]
        optimal = mat_cm["optimal"]
        solverstats = mat_cm["solverstats"]
        exp1 = mat_cm["expS"] # one way to get expansion
        approx1 = exp1/lb
        if optimal
            mm = 50
            eS, expS = round_LP_many(X,mm,Elist)
            approx = expS/lb
            approxbest = min(approx,approx1) # report best of both rounding schemes
            solvetime = solverstats[3]

            push!(approx_lp,approxbest)
            push!(runs_lp,solvetime)

        else
            push!(approx_lp,-1)
            push!(runs_lp,-1)
        end
    catch
        push!(approx_lp,-1)
        push!(runs_lp,-1)
    end

    # Load results from SDP relaxation
    try
        mat_sdp = matread("../data/main-data/mathoverflow/math_$(SetInds[num])_SDP_solution.mat")
        X = mat_sdp["X"]
        X = X+ X'
        totaltime = mat_sdp["timeSDP"]
        setuptime = totaltime - times[num]
        lb = mat_sdp["lowerbound"]
        if solved[num]
            eS, expS = round_hyperLR(X,Elist)
            approx = expS/lb
            push!(approx_sdp,approx)
            push!(runs_sdp,times[num])

        else
            push!(approx_sdp,-1)
            push!(runs_sdp,-1)
        end
    catch
        push!(approx_sdp,-1)
        push!(runs_sdp,-1)
    end

end

approx_cm = [approx_cm1 approx_cm2 approx_cm3]
runs_cm = [runs_cm1 runs_cm2 runs_cm3]
Hypstats = [N M Mu]

matwrite("plotdata_mathoverflow.mat", Dict("Hypstats"=>Hypstats,"approx_cm"=>approx_cm, "runs_cm"=>runs_cm, "approx_lp"=>approx_lp, "runs_lp"=>runs_lp, "approx_sdp"=>approx_sdp, "runs_sdp"=>runs_sdp))