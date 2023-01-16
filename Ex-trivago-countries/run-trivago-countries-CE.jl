using MAT
using SparseArrays
using MatrixNetworks

## Code
include("../src/card-based-CE.jl")

## Run the smart clique expansion method

M = matread("../data/trivago-countries/trivago_countries_large_summary.mat")
Labels = M["SpecialLabels"]
LabelNames = M["LabelNames"]

inds = [5, 8, 17, 20]

for ii = 1:4
i = inds[ii]
countryind = Labels[i]
countryname = LabelNames[countryind]
l = countryind
mat = matread("../data/trivago-countries/trivago_countries_$(l)_2core.mat")
H = mat["H"] 
m,n = size(H)
order =  vec(sum(H,dims = 2))
order = round.(Int64,order)
d = vec(sum(H,dims = 1))
@assert(length(d) == n)
@assert(length(order) == m)
mu = sum(order)
maxr = round(Int64,maximum(order))
meanr = round(sum(order)/length(order),digits = 1)
println("$countryname: \t  $n \t  $m \t  $maxr \t $meanr \t $mu")
volA = sum(d)
Edges = incidence2elist(H);

## Run smart ce
deltas = [1 1.1 1.5 2 5 10 100]
for a = 1:length(deltas)
    delta = deltas[a]
    EdgesW = Delta_EdgesW(order,delta)     # hyperedge splitting functions
    nodeweights = d 

    smartsweep = true
    graphnormalize = true
    S, condS, reducetime, sweeptime, y = Smart_CE_Spectral(Edges,EdgesW,nodeweights,smartsweep,graphnormalize)

    matwrite("Output/CE_tric_$(l)_delta_$(delta)_gnorm_$graphnormalize.mat",Dict("S"=>S,"condS"=>condS,"reducetime"=>reducetime,"sweeptime"=>sweeptime,"y"=>y))
    lS = length(S)
    if lS > n/2
        lS = n-length(S)
    end
    println("$delta \t $lS $condS $(sweeptime+reducetime)")
end
end