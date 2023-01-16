# In this file I will load the mathoverflow dataset and extract a number of small hypergraphs
# (based on Q&A topic) to use as benchmarks for the experiments

using SparseArrays
using LinearAlgebra
include("../../src/hypergraph-helper-functions.jl")

dataset = "mathoverflow-answers"

## Read in the node labels
f = readlines("../raw-data/$dataset/node-labels-$dataset.txt")
n = length(f)
Labels = Vector{Vector{Int64}}()
for i = 1:n
    labels = parse.(Int64,split(f[i],","))
    push!(Labels,labels)
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
f = readlines("../raw-data/$dataset/label-names-$dataset.txt")
t = length(f)
LabelNames = Vector{String}()
for i = 1:t
    push!(LabelNames,f[i])
end

## Store a label matrix, from nodes to labels
Nodes2Labels = elist2incidence(Labels,t)


matwrite("../larger-hypergraphs/mathoverflow-answers-all.mat", Dict("H"=>H,"LabelNames"=>LabelNames,"Nodes2Labels"=>Nodes2Labels))



## Now separate into small pieces
Names = LabelNames
L = Nodes2Labels
SpecialLabels = Vector{Int64}()
Sizes = Vector{Int64}()
for l = 1:t
    qs = findall(x->x==1,L[:,l]) # get all questions on topic l
    num = length(qs)
    if num > 100 && num < 300
        Hq,  dq, orderq = largest_cc_star(H[:,qs])
        mu = sum(orderq)
        m,n = size(Hq)
        maxr = round(Int64,maximum(orderq))
        meanr = round(sum(orderq)/length(orderq),digits = 1)
        push!(SpecialLabels,l)
        push!(Sizes,n)
        if n > 75
            println("$l \t $maxr \t $meanr \t $n \t $m \t $mu \t $(Names[l])")
        end 
    end
end


## Now order them
p = sortperm(Sizes)
reordered = SpecialLabels[p]

for i = 1:length(p)
    l = reordered[i]
    qs = findall(x->x==1,L[:,l]) # get all questions on topic l
    num = length(qs)
    if num > 100 && num < 300
        Hq,  dq, orderq = largest_cc_star(H[:,qs])
        mu = sum(orderq)
        m,n = size(Hq)
        maxr = round(Int64,maximum(orderq))
        meanr = round(sum(orderq)/length(orderq),digits = 1)
        if n > 75
            println("$l \t $maxr \t $meanr \t $n \t $m \t $mu \t $(Names[l])")
        end 
    end
end


## Extract and save a subset of small hypergraphs at various sizes
SetInds = [717; 313; 386; 88;53;147;79;227;140;851;420;690;537;285;362;493;523;180;395;293;89;465;397]
stats = zeros(length(SetInds),5)
topics = Vector{String}()
for i = 1:length(SetInds)
    l = SetInds[i]
    qs = findall(x->x==1,L[:,l]) # get all questions on topic l
    Hq,  dq, orderq = largest_cc_star(H[:,qs])
    mu = sum(orderq)
    m,n = size(Hq)
    maxr = round(Int64,maximum(orderq))
    meanr = round(sum(orderq)/length(orderq),digits = 1)
    stats[i,:] = [n,m,mu,maxr,meanr]
    println("$l \t $maxr \t $meanr \t $n \t $m \t $mu \t $(Names[l])")
    topic = Names[l]
    push!(topics,topic)
    matwrite("../mathoverflow-small/math_$(l)_lcc.mat", Dict("H"=>Hq,"topic"=>topic))

end

## Save a summary of the datasets
matwrite("../mathoverflow-small/math_summary.mat",Dict("SetInds"=>SetInds,"stats"=>stats,"topics"=>topics,"statkey" => "num nodes, num edges, sum of edge sizes, max hyperedge size, mean hyperedge size"))
