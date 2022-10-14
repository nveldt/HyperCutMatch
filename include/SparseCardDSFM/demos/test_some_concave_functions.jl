include("../src/pwl_approx.jl")
tol = 1e-10
k = 20
r = floor(Int64,(k)/2)
Is = collect(0:k)
# w = min.(Is,4)            # delta-linear

# Playing around with numbers, we see that the asymptotically
# bad behavior of sqrt doesn't set in until extremely large hypereges
w = Is .* (k .- Is)       # clique
w = sqrt.(Is)
# w = log2.(Is .+ 1)
w = Is.^(.8)
w = w[1:r+1]
Is = Is[1:r+1]

pl = scatter(Is,w,legend = false)

# next

epsi = 0.0
Lines, M, X, Y, Ranges = get_pwl_approx(w,epsi,false)
xs = 0:0.01:Is[end]
mlast = M[end]

if mlast < 0
    # negative sloped lines not needed
    # replace with flat line
    M[end] = 0
    X[end] = Is[end]
    Y[end] = w[end]
    pop!(Lines)
elseif Ranges[end][1] == r-1 && Ranges[end-1][2] == r-1 && mlast > 0
    # Covered point r-1 twice. Get rid of it
    M[end] = 0
    X[end] = Is[end]
    Y[end] = w[end]
    pop!(Lines)
elseif mlast > 0
    push!(M,0)
    push!(X,Is[end])
    push!(Y,w[end])
end
push!(Lines, x->w[end]*ones(length(x)))

# construct gadget
b = zeros(length(X)-1)
for t = 1:length(X)-1
    x,y = line_intersection(M[t],M[t+1],X[t],X[t+1],Y[t],Y[t+1])
    b[t] = x
end

a = zeros(length(M)-1)
for i = 1:length(a)
    a[i] = M[i]- M[i+1]
end
J = length(a)

p = x -> sum(a[j].*min.(x,b[j]) for j = 1:J)

wapprox = p(0:r)
@show length(Lines)

# plotting
for l in Lines[1:end]
    plot!(pl,xs,l(xs))
end
plot!(pl,ylimit = [0,w[end]+1])
# xs = 0:0.01:r
# ys = p(xs)
# plot!(xs,ys,linewidth = 2)
