using Plots

include("src/pwl_approx.jl")

## Define a submodular cardinality based function via Combined gadget
k = 2
w = random_scb_function(k::Int64)

w = 10*w/maximum(w)
pl = scatter(0:k,w,legend = false, xticks = 0:k)
epsi = 10.0
z0, zk, a, b, cgf = SubCardFun_to_CGF_weights(w,epsi,false)
Jtrue = length(a) + 1
J = Jtrue
xs = 0:0.01:k
ys = cgf(xs)
plot!(pl,xs,ys,color= :blue)

# Refine to get better approximation
z0l, zkl, al, bl, cgfl, best_eps = Refined_SCB_to_CGF(w,epsi)
xs = 0:0.01:k
ys = cgfl(xs)
plot!(pl,xs,ys,color= :red)
