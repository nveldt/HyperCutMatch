using MAT
using Plots
using LaTeXStrings
function best_alphas(alphas)
    n = length(alphas)
    bests = zeros(n)
    for i = 1:n
        bests[i] = minimum(alphas[1:i])
    end
    return bests
end

## Code
include("../src/HyperCutMatching.jl");

## Load dataset
datasets = ["Cities_H", "Zoo_H", "mathoverflow-answers-all","tripadvisor","Amazon9_lcc","trivago-cities-all"]
jj = 5
i = 6
hyper = datasets[i]

M = matread("Output/HCM-$(hyper)-$(jj).mat")
LBs = M["LBs"]
Lams = M["Lams"]
Runtimes = M["RuntimesHCM"]
Approx = M["Approx"]

Alphas = best_alphas(M["Alphas"])
## Plot approximations
lw = 1.5
ms = 6
gfs = 12
tfs = 10
titlesize = 14
s1 = 300
s2 = 200
color0 = RGB(27/255,158/255,119/255)
color4 = RGB(217/255,95/255,2/255)
color1 = :brown
color2 = :orange
color3 = :blue
ms1 = :star2
ms2 = :diamond
ms3 = :circle
title = ""
leg = :best
xlab = "Iteration"
ylab = "Approximation Ratio"
p = plot()
p = Plots.plot(fontfamily = "helvetica",linewidth = lw,yscale = :identity,legend = leg,grid = false, size = (s1,s2),color = :gray)
title!(p,title,titlefontsize = titlesize)
xaxis!(p,xlab,tickfontsize=tfs,guidefontsize = gfs)
yaxis!(p,ylab,tickfontsize=tfs,guidefontsize = gfs)
xs = 20:length(Approx)
ys = Approx[xs]
plot!(p,xs, ys, fillalpha=0.3, color = color1, linewidth = lw,markersize = 0,markershape = ms1, markerstrokecolor = color1, label = "")
savefig("Figures/Approx_$(hyper).pdf")

## Lower bound
leg = :best
xlab = "Iteration"
ylab = "Lower Bound"
p = plot()
p = Plots.plot(fontfamily = "helvetica",linewidth = lw,yscale = :identity,legend = leg,grid = false, size = (s1,s2),color = :gray)
title!(p,title,titlefontsize = titlesize)
xaxis!(p,xlab,tickfontsize=tfs,guidefontsize = gfs)
yaxis!(p,ylab,tickfontsize=tfs,guidefontsize = gfs)
str = 10
xs = str:length(LBs)
ys = LBs[xs]
plot!(p,xs, ys, color = :green, linewidth = lw,markersize = 0,markershape = ms1, markerstrokecolor = color1, label = "")

savefig("Figures/LB_$(hyper).pdf")

## Upper bound
leg = :best
xlab = "Iteration"
ylab = "Best Conductance "
p = plot()
p = Plots.plot(fontfamily = "helvetica",linewidth = lw,yscale = :identity,legend = leg,grid = false, size = (s1,s2),color = :gray)
title!(p,title,titlefontsize = titlesize)
xaxis!(p,xlab,tickfontsize=tfs,guidefontsize = gfs)
yaxis!(p,ylab,tickfontsize=tfs,guidefontsize = gfs)
str = 1
xs = str:length(LBs)
ys = Alphas[xs]
plot!(p,xs, ys, color = :blue, linewidth = lw,legend = false)

savefig("Figures/UB_$(hyper).pdf")

end