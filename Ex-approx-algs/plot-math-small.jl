using MAT
using StatsBase
using Plots
using LaTeXStrings

dataset = "mathoverflow"
M = matread("plotdata_$dataset.mat")
Hyperstats = M["Hypstats"]
approx_cm=M["approx_cm"]
runs_cm=M["runs_cm"]
exp_cm=M["exp_cm"]
approx_lp=M["approx_lp"]
runs_lp=M["runs_lp"]
exp_lp = M["exp_lp"]
approx_sdp=M["approx_sdp"]
runs_sdp=M["runs_sdp"]
exp_sdp = M["exp_sdp"]
bestapp_sdp = M["bestapp_sdp"]
N = Hyperstats[:,1]
M = Hyperstats[:,2]
Mu = Hyperstats[:,3]

keep_alp = findall(x->x>0,approx_lp)
keep_sdp = findall(x->x>0,approx_sdp)


e_lp = exp_lp[keep_alp]
a_lp = approx_lp[keep_alp]
r_lp = runs_lp[keep_alp]
Mu_lp = Mu[keep_alp]
N_lp = N[keep_alp]

e_sdp = exp_sdp[keep_sdp]
ba_sdp = bestapp_sdp[keep_sdp]
a_sdp = approx_sdp[keep_sdp]
Mu_sdp = Mu[keep_sdp]
N_sdp = N[keep_sdp]

r_sdp = runs_sdp[keep_sdp]

## Color scheme
c1 = :blue
c2 = :lightblue
c3 = :orange 
c4 = :red

## Approx vs. N
s1 = 300
s2 = 225
f = plot(title = "",grid = false,legend = :false, xlabel = "n", ylabel = "Approximation Ratio")
scatter!(f,N,approx_cm[:,1],markerstrokewidth = 0, label = "HCM-1", size = (s1,s2),color = c1)
scatter!(f,N,approx_cm[:,2],markerstrokewidth = 0, label = "HCM-2", size = (s1,s2),color = c2)

scatter!(f,N_lp,a_lp,markerstrokewidth = 0, label = "LP",color = c3)
scatter!(f,N_sdp,a_sdp,markerstrokewidth = 0, label = "SDP",color = c4)
savefig("Figures/math-small-approx.pdf")


## Expansion vs. N
s1 = 300
s2 = 225
f = plot(title = "",grid = false,legend = :false, xlabel = "n", ylabel = "Expansion")
scatter!(f,N,exp_cm[:,1],markerstrokewidth = 0, label = "HCM-1", size = (s1,s2),color = c1)
scatter!(f,N,exp_cm[:,2],markerstrokewidth = 0, label = "HCM-2", size = (s1,s2),color = c2)

scatter!(f,N_lp,e_lp,markerstrokewidth = 0, label = "LP",color = c3)
scatter!(f,N_sdp,e_sdp,markerstrokewidth = 0, label = "SDP",color = c4)





## Run vs. N
s1 = 300
s2 = 225
f = plot(title = "",grid = false,yaxis = :log10,legend = false, xlabel = "n", ylabel = "Runtime (s)")
scatter!(f,N,runs_cm[:,1],markerstrokewidth = 0, label = "HCM-1", size = (s1,s2),color =c1)
scatter!(f,N,runs_cm[:,2],markerstrokewidth = 0, label = "HCM-2", size = (s1,s2),color = c2)
scatter!(f,N_lp,r_lp,markerstrokewidth = 0, label = "LP",color = c3)
scatter!(f,N_sdp,r_sdp,markerstrokewidth = 0, label = "SDP",color = c4)

savefig("Figures/math-small-runtimes.pdf")


## Approx vs. Problem Size
s1 = 300
s2 = 225
xlab = L"\sum_{e \in E} |e|"
f = plot(title = "",grid = false,legend = :false, xlabel = xlab, ylabel = "Approx Ratio")
scatter!(f,Mu,approx_cm[:,1],markerstrokewidth = 0, label = "HCM-1", size = (s1,s2),color = c1)
scatter!(f,Mu,approx_cm[:,2],markerstrokewidth = 0, label = "HCM-2", size = (s1,s2),color = c2)
scatter!(f,Mu_lp,a_lp,markerstrokewidth = 0, label = "LP",color = c3)
scatter!(f,Mu_sdp,a_sdp,markerstrokewidth = 0, label = "SDP",color = c4)
# savefig("Figures/math-small-approx-mu.pdf")


## Run vs. Problem Size
s1 = 300
s2 = 225
xlab = L"\sum_{e \in E} |e|"
f = plot(title = "",grid = false,yaxis = :log10,legend = :false, xlabel = xlab, ylabel = "Runtime (s)")
scatter!(f,Mu,runs_cm[:,1],markerstrokewidth = 0, label = "HCM-1", size = (s1,s2),color = c1)
scatter!(f,Mu,runs_cm[:,2],markerstrokewidth = 0, label = "HCM-2", size = (s1,s2),color = c2)
scatter!(f,Mu_lp,r_lp,markerstrokewidth = 0, label = "LP",color = c3)
scatter!(f,Mu_sdp,r_sdp,markerstrokewidth = 0, label = "SDP",color = c4)
# savefig("Figures/math-small-runtimes-mu.pdf")
