using MAT
using SparseArrays
using MatrixNetworks
include("../../src/hypergraph-helper-functions.jl")
mat = matread("../amazon-9/Amazon_9.mat")
H = mat["H"]
labels = mat["labels"]
names = mat["names"]
m,n = size(H)

## This code simply re-orders some labels from 1-9, rather than using original labels from
# a larger hypergraph (Veldt et al. 2020, KDD)
counts = zeros(9)
labelold = unique(labels)
for i = 1:9
    counts[i] = length(findall(x->x==labelold[i],labels))
end
Label2Count = [labelold counts]
p = sortperm(counts)
Label2Name = [Label2Count[p,:] names]
old2new = Dict()
for j = 1:9
    old2new[Label2Name[j,1]] = j
end
nodelabels = zeros(Int64,n)
for i = 1:n
    oldlabel = labels[i]
    newlabel = old2new[labels[i]]
    nodelabels[i] = newlabel
end

## Take largest connected component
A = [spzeros(n,n) sparse(H'); H spzeros(m,m)]
Acc,p = largest_component(A)
pnodes = p[1:n]
pedges = p[n+1:n+m]
H = H[pedges,pnodes]
d = vec(sum(H,dims = 1))
order = vec(sum(H,dims = 2))
keep_edges = findall(x->x>1,order)
H = H[keep_edges,:]
d = vec(sum(H,dims = 1))
order = vec(sum(H,dims = 2))
labels = nodelabels[pnodes]

matwrite("../larger-hypergraphs/Amazon9_lcc.mat", Dict("H"=>H,"labels"=>labels,"names"=>names))


## Extract each individual hypergraph 
mat = matread("../larger-hypergraphs/Amazon9_lcc.mat")
Hall = mat["H"]
labels = mat["labels"]
names = mat["names"]
m,n = size(H)
for i = 1:9
    nodes = findall(x->x==i,labels)
    H = Hall[:,nodes]
    H, d, order = largest_cc_star(H)
    mu = sum(order)
    m,n = size(H)
    maxr = round(Int64,maximum(order))
    meanr = round(sum(order)/length(order),digits = 1)
    println("$maxr \t $meanr \t $n \t $m \t $mu \t $(names[i])")
    matwrite("../amazon-9/Am_$(i)_lcc.mat", Dict("H"=>H))
end
