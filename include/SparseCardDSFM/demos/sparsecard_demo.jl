using SparseArrays

include("../hypergraphs.jl")
include("../SparseCard.jl")
include("../pwl_approx.jl")

## The HyperModualrity package includes a number of hypergraph datasets we can use for testing
using HyperModularity 
hypermodularity_datasets() 

## Load dataset and solve
dataset = "contact-primary-school-classes"
# dataset = "mathoverflow-answers"
maxsize = 25	# max hyperedge size
minsize = 2	# min hyperedge size
return_labels = true
Edges, Labels = read_hypergraph_data(dataset,maxsize,minsize,return_labels)
pb = 0.2:0.1:0.3
EdgesW = RandomSplittingFunctions(Edges,pb)
P = sortperm(Labels)
n = length(Labels)
refined = false
epsilon = 0.1*ones(maxsize)


## Solve the problem with two different max-flow subroutines
eS_pr, cutval_pr, approx_pr, flowtime_pr, reduce_time, worst_eps = SparseCard(Edges,EdgesW,n,epsilon;refined = false)
eS_gur, cutval_gur, approx_gur,flowtime_gur, reduce_time, worst_eps = SparseCard(Edges,EdgesW,n,epsilon;flowsolver = "gurobi")


## print output details

println("Reduce time = $reduce_time \t Approx =  $worst_eps")
println("PushRelabel: $cutval_pr \t  $approx_pr \t $flowtime_pr \t")
println("Gurobi-Flow: $cutval_pr \t  $approx_pr \t $flowtime_gur")

## More detail for reduction
A, svec, tvec = CardDSFM_reduction(Edges,EdgesW,n,epsilon,refined)
N = size(A,1)

## Solve the max-flow with push-relabel code
F = maxflow(A,svec,tvec,0)
Sall = F.source_nodes[2:end].-1
ct = F.cutvalue
Small = source_nodes_min(F,1e-10)[2:end].-1

S = intersect(Sall,1:n)
eS = zeros(n)
eS[S] .= 1
cutS = eval_fS(Edges,EdgesW,eS)

## Find the best cut with certain nodes fixed
T = setdiff(1:n,S)
S2, cutval = eval_min_dircut(A,svec,tvec,S,T)

## Solve the max-flow with Gurobi
Svec, flowval, f = Gubori_MaxFlow(A,svec,tvec,1)
eS = Svec[1:n]
cut_gurobi = eval_fS(Edges,EdgesW,Svec)
