using SparseArrays
using LinearAlgebra

include("../../src/hypergraph-helper-functions.jl")
dataset = "TrivagoClickout"

## Read in the node labels
f = readlines("../raw-data/$dataset/node-labels-$dataset.txt")
n = length(f)
Labels = Vector{Int64}()
for i = 1:n
    push!(Labels, parse(Int64, f[i]))
end

## Read in the hypergraph 
f = readlines("../raw-data/$dataset/hyperedges-$dataset.txt")
m = length(f)
Edges = Vector{Vector{Int64}}()
for i = 1:m
    edge = parse.(Int64,split(f[i],","))
    push!(Edges,edge)
end
H = elist2incidence(Edges,n)

## Read in the label names
f   = readlines("../raw-data/$dataset/label-names-$dataset.txt")
t = length(f)
LabelNames = Vector{String}()
for i = 1:t
    push!(LabelNames,f[i])
end
 
## Extract larger hypergraphs
SpecialLabels = Vector{Int64}()
for l = 1:t
    qs = findall(x->x==l,Labels) # get accomodations with Label t
    num = length(qs)
    if  num > 1000
        Hq,  dq, orderq = largest_cc_star(H[:,qs])
        mu = sum(orderq)
        m,n = size(Hq)
        maxr = round(Int64,maximum(orderq))
        meanr = round(sum(orderq)/length(orderq),digits = 1)
        push!(SpecialLabels,l)
        Hcore, d, order = hyperkcore_multi(Hq,2)
        if in(l,[6;10;23;26])
            # Just saving the 4 datasets were are interested in
            matwrite("../trivago-countries/trivago_countries_$(l)_2core.mat", Dict("H"=>Hcore))
        end
        if n > 75
            println("$l \t $maxr \t $meanr \t $n \t $m \t $mu \t $(LabelNames[l])")
        end 
    end
end

## Save the data
matwrite("../trivago-countries/trivago_countries_large_summary.mat",Dict("SpecialLabels"=>SpecialLabels,"LabelNames"=>LabelNames))


