using SparseArrays
using LinearAlgebra
include("../../src/hypergraph-helper-functions.jl")

dataset = "TrivagoClickout"

cities = true
## Read in the node labels
if cities
    f = readlines("../raw-data/$dataset/node-labels-$dataset-cities.txt")
else
    f = readlines("../raw-data/$dataset/node-labels-$dataset.txt")
end
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
if cities
    f   = readlines("../raw-data/$dataset/label-names-$dataset-cities.txt")
else
    f   = readlines("../raw-data/$dataset/label-names-$dataset.txt")
end
t = length(f)
LabelNames = Vector{String}()
for i = 1:t
    push!(LabelNames,f[i])
end

if cities
    matwrite("../larger-hypergraphs/trivago-cities-all.mat", Dict("H"=>H,"LabelNames"=>LabelNames,"Labels"=>Labels))
else
    matwrite("../larger-hypergraphs/trivago-countries-all.mat", Dict("H"=>H,"LabelNames"=>LabelNames,"Labels"=>Labels))
end 

## Now separate into small pieces
SpecialLabels = Vector{Int64}()
Sizes = Vector{Int64}()
for l = 1:t
    qs = findall(x->x==l,Labels) # get accomodations with Label t
    num = length(qs)
    if num < 300 && num > 100
        Hq,  dq, orderq = largest_cc_star(H[:,qs])
        mu = sum(orderq)
        m,n = size(Hq)
        maxr = round(Int64,maximum(orderq))
        meanr = round(sum(orderq)/length(orderq),digits = 1)
        push!(SpecialLabels,l)
        push!(Sizes,n)
        if n > 75
            println("$l \t $maxr \t $meanr \t $n \t $m \t $mu \t $(LabelNames[l])")
        end 
    end
end


## Now order them
p = sortperm(Sizes)
reordered = SpecialLabels[p]

for i = 1:length(p)
    l = reordered[i]
    qs = findall(x->x==l,Labels) # get accomodations with Label t
    num = length(qs)
    if num > 100 && num < 300
        Hq,  dq, orderq = largest_cc_star(H[:,qs])
        mu = sum(orderq)
        m,n = size(Hq)
        maxr = round(Int64,maximum(orderq))
        meanr = round(sum(orderq)/length(orderq),digits = 1)
        if n > 75
            println("$l \t $maxr \t $meanr \t $n \t $m \t $mu \t $(LabelNames[l])")
        end 
    end
end


## Extract a subset of small hypergraphs, with varying sizes, to use as benchmarks

SetInds = [257, 1007, 530, 50, 170, 515, 549, 8, 2160, 2718, 742, 1589, 1019, 993, 444, 1730, 961, 341, 299, 1365, 647, 667, 651, 765, 525, 121, 1225, 442, 284, 1014, 155, 518, 119, 88, 614, 860, 1030, 125, 559, 379, 105]
stats = zeros(length(SetInds),5)
cities = Vector{String}()
for i = 1:length(SetInds)
    l = SetInds[i]
    qs = findall(x->x==l,Labels) # get all questions on topic l
    Hq,  dq, orderq = largest_cc_star(H[:,qs])
    mu = sum(orderq)
    m,n = size(Hq)
    maxr = round(Int64,maximum(orderq))
    meanr = round(sum(orderq)/length(orderq),digits = 1)
    stats[i,:] = [n,m,mu,maxr,meanr]
    println("$l \t $maxr \t $meanr \t $n \t $m \t $mu \t $(LabelNames[l])")
    city = LabelNames[l]
    push!(cities,city)
    matwrite("../trivago-small/trivago_$(l)_lcc.mat", Dict("H"=>Hq,"city"=>city))

end

## Save a summary of the stored hypergraphs

matwrite("../trivago-small/trivago_summary.mat",Dict("SetInds"=>SetInds,"stats"=>stats,"cities"=>cities,"statkey" => "num nodes, num edges, sum of edge sizes, max hyperedge size, mean hyperedge size"))
