using MAT
using Plots
using LaTeXStrings
datasets = ["Newsgroups", "Mushrooms", "Covertype45", "Covertype67"]

i = 4
name = datasets[i]
gnormstring = true
M = matread("plotdata_benchmark_gnorm_true.mat")
hcm_cm = M["HCM_conds_mean"][:,i]
hcm_cs = M["HCM_conds_std"][:,i]
hcm_tm = M["HCM_time_mean"][:,i]
hcm_ts = M["HCM_time_std"][:,i]
hcm_am = M["HCM_app_mean"][:,i]
hcm_as = M["HCM_app_std"][:,i]
ipm_c = M["IPM_conds"][:,i]
ipm_t = M["IPM_time"][:,i]
ipm_am = M["IPM_app_mean"][:,i]
ipm_as = M["IPM_app_std"][:,i]
ce_c = M["CE_conds"][:,i]
ce_t = M["CE_time"][:,i]
ce_am = M["CE_app_mean"][:,i]
ce_as = M["CE_app_std"][:,i]

alphas = collect(round.(LinRange(0.000,0.04,11),digits = 3))

## Plot approximations
lw = 1.5
ms = 6
gfs = 12
tfs = 10
titlesize = 14
s1 = 300
s2 = 200
sc = 2
color1 = RGB(27/255,158/255,119/255)
color2 = RGB(217/255,95/255,2/255)
color1 = :brown
color2 = :orange
color3 = :blue
ms1 = :star2
ms2 = :diamond
ms3 = :circle
title = ""
leg = false
if i == 4
    leg = :bottomright
end
xlab = L"\alpha"
ylab = "Approximation Ratio"
ylims = (minimum(hcm_am)-.1,maximum(ipm_am)+.1)
p = plot()
p = Plots.plot(fontfamily = "helvetica",linewidth = lw,yscale = :identity,legend = leg,grid = false, size = (s1,s2),color = :gray)
title!(p,title,titlefontsize = titlesize, xticks = (1:2:11,["0", "0.008", "0.016", "0.024", "0.032", "0.4"])) #,ylim = ylims)
xaxis!(p,xlab,tickfontsize=tfs,guidefontsize = gfs)
yaxis!(p,ylab,tickfontsize=tfs,guidefontsize = gfs)

plot!(p,1:11,[ce_am ce_am], fillrange=[ce_am-ce_as ce_am+ce_as], fillalpha=0.3, color = color1, linewidth = lw, markersize = ms,
    markershape = ms1, markerstrokecolor = color1, label = "")
plot!(p,1:11,[ipm_am ipm_am], fillrange=[ipm_am-ipm_as ipm_am+ipm_as], fillalpha=0.3, color = color2, linewidth = lw, markersize = ms,
    markershape = ms2, markerstrokecolor = color2, label = "")
plot!(p,1:11,[hcm_am hcm_am], fillrange=[hcm_am-hcm_as hcm_am+hcm_as], fillalpha=0.3, color = color3, linewidth = lw, markersize = ms,
    markershape = ms3, markerstrokecolor = color3, label = "")

plot!(p,1:11,ce_am, color = color1, linewidth = lw, markersize = ms,
    markershape = ms1, markerstrokecolor = color1, label = "CE+HCM")
plot!(p,1:11,ipm_am, fillalpha=0.3, color = color2, linewidth = lw, markersize = ms,
    markershape = ms2, markerstrokecolor = color2, label = "IPM+HCM")
plot!(p,1:11,hcm_am, fillalpha=0.3, color = color3, linewidth = lw, markersize = ms,
    markershape = ms3, markerstrokecolor = color3, label = "HCM")

# savefig("Figures/Ratios_$(name)_$gnormstring.pdf")

##
leg = false
ylab = "Runtime (s)"
p = plot()
p = Plots.plot(fontfamily = "helvetica",linewidth = 2,yscale = :log10,legend = leg,grid = false, size = (s1,s2),color = :gray)
title!(p,title,titlefontsize = titlesize,xticks = (1:2:11,["0", "0.008", "0.016", "0.024", "0.032", "0.4"]))
xaxis!(p,xlab,tickfontsize=tfs,guidefontsize = gfs)
yaxis!(p,ylab,tickfontsize=tfs,guidefontsize = gfs)

plot!(p,1:11,ce_t, color = color1, linewidth = lw, markersize = ms,
    markershape = ms1, markerstrokecolor = color1, label = "")
plot!(p,1:11,ipm_t, fillalpha=0.3, color = color2, linewidth = lw, markersize = ms,
    markershape = ms2, markerstrokecolor = color2, label = "")
plot!(p,1:11,[hcm_tm hcm_tm], fillrange=[hcm_tm-hcm_ts hcm_tm+hcm_ts], fillalpha=0.3, color = color3, linewidth = lw, markersize = ms,
    markershape = ms3, markerstrokecolor = color3, label = "")

savefig("Figures/Runtimes_$(name)_$gnormstring.pdf")


