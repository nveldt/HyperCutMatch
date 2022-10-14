using MAT
using StatsBase
using Statistics

## Load information about the four datasets
M = matread("../data/trivago-large/trivago_countries_large_summary.mat")
Labels = M["SpecialLabels"]
LabelNames = M["LabelNames"]

inds = [5, 8, 17, 20]
orig_labels = Labels[inds]
Names = LabelNames[orig_labels]

## Load and plot
deltas = [1 1.1 1.5 2 5 10 100]
dts = length(deltas)
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



for i = 1:4

l = orig_labels[i]
countryname = LabelNames[l]
mat = matread("../data/trivago-large/trivago_countries_$(l)_2core.mat")
H = mat["H"] 
order = round.(Int64,vec(sum(H,dims = 2)))
d = vec(sum(H,dims = 1))
volA = sum(d)
m,n = size(H)
Edges = incidence2elist(H);

numtimes = 5
hcm_times = zeros(length(deltas),numtimes)
hcm_conds = zeros(length(deltas),numtimes)
Lbs = zeros(length(deltas),numtimes)
hcm_approx = zeros(length(deltas),numtimes)
ipm_times = zeros(length(deltas))
ipm_conds = zeros(length(deltas))
ipm_approx = zeros(length(deltas),numtimes)
ce_times = zeros(length(deltas))
ce_conds = zeros(length(deltas))
ce_approx = zeros(length(deltas),numtimes)


for a = 1:length(deltas)
    delta = deltas[a]
    EdgesW = Delta_EdgesW(order,delta)    

    ## Load IPM output
    if isinteger(delta)
        del = round(Int64,delta)
        m2 = matread("Output/IPM_tric_$(l)_2core_$(del)_1_combinewarms.mat")
    else
        m2 = matread("Output/IPM_tric_$(l)_2core_$(delta)_1_combinewarms.mat")
    end
    eipm = vec(m2["eipm"])
    Sipm = findall(x->x>0,eipm)
    ipmtime = m2["ipmtime"]
    ipmcond = m2["ipmcond"]
    condS, volS, cutS = tl_cond(H,Sipm,d,delta,volA,order)
    @assert(abs(ipmcond-condS) < 1e-12)
    ipm_times[a] = ipmtime
    ipm_conds[a] = ipmcond

    ## Load CE data
    if isinteger(delta)
        del = round(Int64,delta)
        m2 = matread("Output/CE_tric_$(l)_2core_$(del)_1_degreewarm.mat")
    else
        m2 = matread("Output/CE_tric_$(l)_2core_$(delta)_1_degreewarm.mat")
    end
    ece = vec(m2["ece"])
    Sce = findall(x->x>0,ece)
    cetime = m2["cetime"]
    cecond = m2["ceCond"]
    condS, volS, cutS = tl_cond(H,Sce,d,delta,volA,order)
    @assert(abs(cecond-condS) < 1e-12)
    ce_times[a] = cetime
    ce_conds[a] = cecond

    ## Load HCM
    for jj = 1:numtimes
        m1 = matread("Output/HCM_tric_$(l)_2core_$(delta)_$(jj)_longer.mat")
        Runs = m1["RuntimesHCM"]
        hcmtime = m1["solvetime"]
        S = m1["S"]
        hcmcond_check = minimum(m1["Conds"])
        hcmcond, volS1, cutS1 = tl_cond(H,S,d,delta,volA,order)
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

matwrite("plotdata_trivago_countries_combinewarms.mat", Dict(
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