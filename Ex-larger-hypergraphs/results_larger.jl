using MAT
using StatsBase
using Statistics

datasets = ["mathoverflow-answers-all","tripadvisor","Amazon9_lcc"]

i = 2
hyper = datasets[i]
mat = matread("../data/larger/$hyper.mat")
H = mat["H"]
H, d, order = largest_cc_star(H)
order = round.(Int64,order)
mu = sum(order)
m,n = size(H)
maxr = round(Int64,maximum(order))
meanr = round(sum(order)/length(order),digits = 1)
# print("$(hyper[1:4]) \t  $n \t  $m \t  $meanr")


runs = zeros(5)
conds = zeros(5)
apps = zeros(5)
lbs = zeros(5)
for jj = 1:5
    M = matread("Output/HCM-$(hyper)-1-$(jj).mat")
    LBs = M["LBs"]
    Lams = M["Lams"]
    Runtimes = M["RuntimesHCM"]
    S = M["S"]
    condS = minimum(M["Alphas"])
    Approx = M["ApproxFac"]
    totaltime = sum(Runtimes)
    runs[jj] = totaltime
    conds[jj] = condS
    apps[jj] = Approx
    lbs[jj] = maximum(LBs)
end

d1 = 3
d2 = 1
cm = round(mean(conds),digits = d1)
cs = round(std(conds),digits = d1)
lm = round(mean(lbs),digits = d1)
ls = round(std(lbs),digits = d1)
rm = round(mean(runs),digits = d2)
rs = round(std(runs),digits = d2)
am = round(mean(apps),digits = 2)
as = round(std(apps),digits = d1)

# println("$(hyper[1:4])&  $n &  $m &  $meanr & $lm {\\small \$\\pm $ls\$} & $cm {\\small \$\\pm $cs\$}  & $am {\\small \$\\pm $as\$ & $rm {\\small \$\\pm $rs\$}}")
println("$(hyper)&  $n &  $m &  $meanr & $am {\\small \$\\pm $as\$ & $rm {\\small \$\\pm $rs\$}}")