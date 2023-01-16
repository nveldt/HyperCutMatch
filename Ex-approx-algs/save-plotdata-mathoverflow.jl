using MAT
using StatsBase
using LaTeXStrings

include("../src/hypergraph-sc-sdp.jl")

M = matread("../data/mathoverflow-small/math_summary.mat")
topics = M["topics"]
stats = M["stats"]
SetInds = M["SetInds"]
numtimes = 10

exp_cm1 = Vector{Float64}()
approx_cm1 = Vector{Float64}()
runs_cm1 = Vector{Float64}()

exp_cm2 = Vector{Float64}()
approx_cm2 = Vector{Float64}()
runs_cm2 = Vector{Float64}()

exp_lp = Vector{Float64}()
approx_lp = Vector{Float64}()
runs_lp = Vector{Float64}()

exp_sdp = Vector{Float64}()
bestapp_sdp = Vector{Float64}()
approx_sdp = Vector{Float64}()
runs_sdp = Vector{Float64}()

N = Vector{Int64}()
M = Vector{Int64}()
Mu = Vector{Int64}()

## Get one for testing

num = 10
mat = matread("../data/mathoverflow-small/math_$(SetInds[num])_lcc.mat")
    H = mat["H"]
    m,n = size(H)
    mu = round(Int64,sum(H))
    Elist = incidence2elist(H)
mat_sdp = matread("Output/mathoverflow/math_$(SetInds[num])_SDP_solution.mat")
X = Symmetric(mat_sdp["X"])
D = decompose_X(X)
eS, expS = round_hyperLR(X,Elist)
eS2, expS2 = round_hyperLR(D,Elist)

##

for num = 1:length(SetInds)
    # Load the data
    mat = matread("../data/mathoverflow-small/math_$(SetInds[num])_lcc.mat")
    H = mat["H"]
    m,n = size(H)
    mu = round(Int64,sum(H))
    Elist = incidence2elist(H)
    push!(N,n)
    push!(M,m)
    push!(Mu,mu)

    # load cut matching results
    mat_cm = matread("Output/mathoverflow/math_$(SetInds[num])_HCM_solution.mat")
    Approx = mat_cm["Approxs"]
    Runtimes = mat_cm["Runtimes"]
    Bounds = mat_cm["Bounds"]
    Exps = Approx.*Bounds
    push!(approx_cm1,mean(Approx[1,:]))
    push!(approx_cm2,mean(Approx[2,:]))
    push!(runs_cm1,mean(Runtimes[1,:]))
    push!(runs_cm2,mean(Runtimes[2,:]))
    push!(exp_cm1, mean(Exps[1,:]))
    push!(exp_cm2, mean(Exps[2,:]))
    @assert(length(Approx[1,:]) == 10)
    minexp_hcm = min(minimum(Exps[1,:]),minimum(Exps[2,:]))

    mat_lp = matread("Output/mathoverflow/mathoverflow_$(SetInds[num])_LP_solution.mat")
    lb = mat_lp["lb"]
    optimal = mat_lp["optimal"]

    if optimal
        solverstats = mat_lp["solverstats"]
        exp1 = mat_lp["expS"] # one way to get expansion
        approx1 = exp1/lb
        solvetime =  solverstats[3]
        push!(approx_lp,approx1)
        push!(runs_lp,solvetime)
        push!(exp_lp,exp1)
    else
        println("LP not optimal")
        push!(approx_lp,-1)
        push!(runs_lp,-1)
        push!(exp_lp,-1)
    end

    # Load results from SDP relaxation
    try
        mat_sdp = matread("Output/mathoverflow/math_$(SetInds[num])_SDP_solution.mat")
        X = mat_sdp["X"]
        totaltime = mat_sdp["timeSDP"]
        lb = mat_sdp["lowerbound"]
        status = mat_sdp["status"]

        if status == "Solved"
            D = decompose_X(X)
            eS1, expS1 = round_hyperLR(X,Elist)
            eS2, expS2 = round_hyperLR(D,Elist)
            expS = min(expS1,expS2)
            approx = expS/lb
            bestexp = minimum([expS; exp1; minexp_hcm])
            bestapprox = bestexp/lb
            push!(approx_sdp,approx)
            push!(runs_sdp,totaltime)
            push!(exp_sdp, expS)
            push!(bestapp_sdp, bestapprox)
        else
            println("SDP finished, but not solved")
            @show SetInds[num]
            push!(approx_sdp,-1)
            push!(runs_sdp,-1)
            push!(exp_sdp,-1)
            push!(bestapp_sdp, -1)
        end

    catch
        push!(approx_sdp,-1)
        push!(runs_sdp,-1)
        push!(exp_sdp,-1)
        push!(bestapp_sdp, -1)
    end
end

exp_cm = [exp_cm1 exp_cm2]
approx_cm = [approx_cm1 approx_cm2]
runs_cm = [runs_cm1 runs_cm2]
Hypstats = [N M Mu]


##
matwrite("plotdata_mathoverflow.mat", Dict("Hypstats"=>Hypstats,
"approx_cm"=>approx_cm, "runs_cm"=>runs_cm, "approx_lp"=>approx_lp, 
"runs_lp"=>runs_lp, "approx_sdp"=>approx_sdp, "runs_sdp"=>runs_sdp,
"exp_cm"=>exp_cm,"exp_sdp"=>exp_sdp,"exp_lp"=>exp_lp,"bestapp_sdp"=>bestapp_sdp))