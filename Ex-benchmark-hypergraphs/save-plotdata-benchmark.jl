using MAT
using StatsBase
using Statistics

include("../src/hypergraph-helper-functions.jl")

gnormstring = "gnorm_true"

## Store output

alphas = collect(round.(LinRange(0.000,0.04,11),digits = 3))

dts = length(alphas)
CE_conds = zeros(dts,4)
CE_times = zeros(dts,4)
CE_app_mean = zeros(dts,4)
CE_app_std = zeros(dts,4)

IPM_conds = zeros(dts,4)
IPM_times = zeros(dts,4)
IPM_app_mean = zeros(dts,4)
IPM_app_std = zeros(dts,4)

HCM_conds_mean = zeros(dts,4)
HCM_conds_std = zeros(dts,4)
HCM_time_mean = zeros(dts,4)
HCM_time_std = zeros(dts,4)
HCM_time_mean = zeros(dts,4)
HCM_app_mean = zeros(dts,4)
HCM_app_std = zeros(dts,4)


## Load and save data

datasets = ["Newsgroups", "Mushrooms", "Covertype45", "Covertype67"]
for i = 1:4
dataset = datasets[i]
mat = matread("../data/benchmark-hypergraphs/$(dataset)_H.mat")
H = mat["H"]
H, d, order = largest_cc_star(H)
order = round.(Int64,order)
mu = sum(order)
m,n = size(H)
maxr = round(Int64,maximum(order))
meanr = round(sum(order)/length(order),digits = 1)
println("\\emph{$(dataset)} &  $n & $m & $meanr & $mu \\\\ ")
volA = sum(d)
Edges = incidence2elist(H)


##
numtimes = 5
hcm_times = zeros(length(alphas),numtimes)
hcm_conds = zeros(length(alphas),numtimes)
Lbs = zeros(length(alphas),numtimes)
hcm_approx = zeros(length(alphas),numtimes)
ipm_times = zeros(length(alphas))
ipm_conds = zeros(length(alphas))
ipm_approx = zeros(length(alphas),numtimes)
ce_times = zeros(length(alphas))
ce_conds = zeros(length(alphas))
ce_approx = zeros(length(alphas),numtimes)

##
for a = 1:length(alphas)
    alpha = alphas[a]
    EdgesW = Alpha_EdgesW(order,alpha)  
    nodeweights = generalized_degree(Edges,EdgesW,n)
    volA = sum(nodeweights)

    ## Load IPM output
    if isinteger(alpha)
        alp = round(Int64,alpha)
        m2 = matread("Output/IPM_$(dataset)_alpha_$(alp)_$(gnormstring).mat")
    else
        m2 = matread("Output/IPM_$(dataset)_alpha_$(alpha)_$(gnormstring).mat")
    end
    eipm = vec(m2["eipm"])
    Sipm = findall(x->x>0,eipm)
    ipmtime = m2["ipmtime"]
    ipmcond = m2["ipmcond"]
    condS = gen_ratio_cut(Edges,EdgesW,Sipm,nodeweights,n,volA)
    @assert(abs(ipmcond-condS) < 1e-12)
    ipm_times[a] = ipmtime
    ipm_conds[a] = ipmcond

    ## Load CE data
    m2 = matread("Output/CE_$(dataset)_alpha_$(alpha)_$(gnormstring).mat")
    Sce = vec(m2["S"])
    reducetime = m2["reducetime"]
    sweeptime = m2["sweeptime"]
    cetime = reducetime+sweeptime
    cecond = m2["condS"]
    condS = gen_ratio_cut(Edges,EdgesW,Sce,nodeweights,n,volA)
    @assert(abs(cecond-condS) < 1e-12)
    ce_times[a] = cetime
    ce_conds[a] = cecond

    ## Load HCM
    for jj = 1:numtimes
        m1 = matread("Output/HCM_$(dataset)_alpha_$(alpha)_$(jj).mat")
        Runs = m1["RuntimesHCM"]
        hcmtime = m1["solvetime"]
        S = m1["S"]
        hcmcond_check = minimum(m1["Conds"])
        hcmcond = gen_ratio_cut(Edges,EdgesW,S,nodeweights,n,volA)

        hcm_times[a,jj] = hcmtime
        hcm_conds[a,jj] = hcmcond
        lb = m1["Lbound"]
        Lbs[a,jj] = lb
        hcm_approx[a,jj] = hcmcond/lb
        ipm_approx[a,jj] = ipmcond/lb
        ce_approx[a,jj] = cecond/lb
        
    end

end


CE_conds[:,i] = ce_conds
CE_times[:,i] = ce_times
CE_app_mean[:,i] = mean(ce_approx,dims = 2)
CE_app_std[:,i] = std(ce_approx,dims = 2)

IPM_conds[:,i] = ipm_conds
IPM_times[:,i] = ipm_times
IPM_app_mean[:,i] = mean(ipm_approx,dims = 2)
IPM_app_std[:,i] = std(ipm_approx,dims = 2)

HCM_conds_mean[:,i] = mean(hcm_conds, dims = 2)
HCM_conds_std[:,i] = std(hcm_conds, dims = 2)
HCM_time_mean[:,i] = mean(hcm_times, dims = 2)
HCM_time_std[:,i] = std(hcm_times, dims = 2)
HCM_app_mean[:,i] = mean(hcm_approx, dims = 2)
HCM_app_std[:,i] = std(hcm_approx, dims = 2)

end

## Save it

matwrite("plotdata_benchmark_$(gnormstring).mat", Dict(
"HCM_conds_mean"=>HCM_conds_mean,
"HCM_conds_std"=>HCM_conds_std,
"HCM_time_mean"=>HCM_time_mean,
"HCM_time_std"=>HCM_time_std,
"HCM_app_mean"=>HCM_app_mean,
"HCM_app_std"=>HCM_app_std,
"IPM_conds"=>IPM_conds,
"IPM_time"=>IPM_times,
"IPM_app_mean"=>IPM_app_mean,
"IPM_app_std"=>IPM_app_std,
"CE_conds"=>CE_conds,
"CE_time"=>CE_times,
"CE_app_mean"=>CE_app_mean,
"CE_app_std"=>CE_app_std))