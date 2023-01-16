# This is an implementation of the inhomogeneous hypergraph clustering algorithm of Li and Milenkovic (NeurIPS 2017),
# specifically for cardinality-based hypergraph splitting functions.
# Previous implementations were in Matlab and focused just on special cases. This implementation is more
# general and faster.

using MAT
using SparseArrays
using LinearAlgebra
using KrylovKit

include("hypergraph-helper-functions.jl")

"""
Implementation of the clique expansion (CE) hypergraph clustering technique of Li and Milenkovic (NeurIPS 2017) 
specifically for cardinality-based hyperedge splitting functions.

Performs smart clique expansion + spectral clustering in order to minimize ratio cuts for
a hypergraph with general hypergraph penalties.

The performance of clique expansion + spectral approaches will depend on several very subtle design choices in
    (1) How the clique expansion is performed (i.e., how to define weights in projection)
    (2) What notion of hypergraph degree one uses, e.g., number of adjaceny hyperedges vs. sum_{e: v in e} w_e({v}) 
        (Furthermore, there are multiple places where degrees show up, so you could use different notions in different places, or not)
    (3) How spectral clustering is performed (sweep cut on the reduced graph vs. sweep cut on original hypergraph) 

This code follows all details as in (Li and Milenkovic, NeurIPS 2017), along with a few added options
to try to improve the hypergraph ratio cut scores when performing the sweep cut.

Inputs: 
    Edges and EdgesW: edgelist and splitting functions
    nodeweights: the hypergraph nodeweight vector to be used for sweepcuts
    
    Different sweep cut versions:
        smartsweep == true: do sweep cut checking generalized ratio cut scores

        graphnormalize == true: normalize the eigenvector based on graph's node degrees
        graphnormalize == false: normalize eigenvector based on hypergraph nodeweights vector

        if smartsweep == false, then graphnormalize = true, since everything is done based only on the reduced graph structure
"""
function Smart_CE_Spectral(Edges::Vector{Vector{Int64}}, EdgesW::Vector{Vector{Float64}},nodeweights::Vector{Float64},smartsweep::Bool = true, graphnormalize = true)
    n = length(nodeweights)    
    volA = sum(nodeweights)
    tic = time()
    A,p = MinDistortionCliqueExpansion(Edges, EdgesW,n)
    reducetime = time()-tic

    tic = time()
    if smartsweep
        # Get second smallest eigenvector of graph normalized Laplacian
        L,Dhalf = nLaplacian(A,true)
        Vl,Vc,convinfo = eigsolve(L + n*LinearAlgebra.I, 2, :SR; tol = 1e-8, maxiter = 1000, verbosity = 0)
        x = Real.(Vc[2])

        # Choose how to degree-normalize the eigenvector
        if graphnormalize
            # use degrees in reduced graph
            y = Dhalf*x
        else
            # use original hypergraph degrees
            Dhalf = Diagonal(nodeweights.^(-1/2))
            y = Dhalf*x
        end

        # "Smart sweep" finds best cut for original hypergraph objective
        S, ratioS = smart_sweep(Edges,EdgesW,nodeweights,y,volA)

    else
        # After reducing to a graph, just do spectral clustering there,
        # ignoring the original hypergraph cut function and nodeweights
        S, condS, y, lam2 = spectral_clustering(A)

        # Compute the hypergraph ratio cut of the set found
        ratioS = gen_ratio_cut(Edges,EdgesW,S,nodeweights,n,volA) 
    end
    sweeptime = time()-tic

    return S, ratioS, reducetime, sweeptime, y

end

"""
Smartsweep: Given a scaled eigenvector of the normalized Laplacian for a the min-distortion reduced graph,
this function performs a sweep cut to return the set that has the best hypergraph ratio score.
"""
function smart_sweep(Edges::Vector{Vector{Int64}}, EdgesW::Vector{Vector{Float64}},nodeweights::Vector{Float64},y::Vector{Float64},volA::Float64)

    n = length(nodeweights)
    p = sortperm(y,rev = true)
    bestcond = 1
    bestS = [p[1]]
    for i = 2:n-1
        S = p[1:i]
        condS = gen_ratio_cut(Edges,EdgesW,S,nodeweights,n,volA)
        if condS < bestcond
            bestS = S
            bestcond = condS
        end
    end
    return bestS, bestcond

end

"""
This performs a minimum distortion clique expansion where the weight of edges in the clique is determined to minimize the 
distortion in the approximation factor. The distortion is the maximum ratio between the projected clique cut penalties
and the hyperege splitting penalty (this is computed locally for each hyperedge).

    Inputs:
        Edges is the set of hyperedges
        EdgesW gives splitting penalties for each hyperedge

    Output:
        A = clique expansion graph
        maxdist = distortion in approximation factor from the edge projection
"""
function MinDistortionCliqueExpansion(Edges::Vector{Vector{Int64}}, EdgesW::Vector{Vector{Float64}},n::Int64)

    I = Vector{Int64}()
    J = Vector{Int64}()
    Vals = Vector{Float64}()
    m = length(Edges)
    maxdist = 1
    for e = 1:m
        edge = Edges[e]
        k = length(edge)

        # Find the projection weight c, which is c = max_i  penalties[i]/ (i*(k-i))
        r = Int64(floor(k/2))
        cliquepen = collect(1:k-1) .* collect(k-1:-1:1)
        cliquepen = cliquepen[1:r]      # clique penalties
        origpen = EdgesW[e][2:end]      # original penalties
        ratio = origpen./cliquepen    
        c = maximum(ratio)              # This guarantees the clique penalties are an upper bound on the original penalties
        distortion = c/minimum(ratio)   # This bounds how much larger the clique penalties can ever be

        if distortion > maxdist
            maxdist = distortion
        end

        for ii = 1:length(edge)
            i = edge[ii]
            for jj = ii+1:length(edge)
                j = edge[jj]
                push!(I,i)
                push!(J,j)
                push!(Vals,c)
            end
        end
    end

    A = sparse(I,J,Vals,n,n)

    return sparse(A+A'),maxdist
end


"""
Return the normalized Laplacian matrix for a graph with adjacency matrix A
"""
function nLaplacian(A::SparseMatrixCSC,dhalf = false)
    d = vec(sum(A,dims = 2))
    Dhalf = Diagonal(d.^(-1/2))
    L = I - Dhalf*A*Dhalf
    if dhalf
        return L, Dhalf
    else
        return L
    end
end

"""
Return the second smallest eigenvalue of the normalized Laplacian matrix 
for the graph whose adjacency matrix is A. 

Simply uses KrylovKit for the eigenvector computation. Could also be approximated in other ways.

Input
    A: n x n adjacency matrix for a connected simple graph

Output

    lam2: second smallest eigenvalue 
    y: vector satisfying lam2 = y'*L*y/(y'*D*y) where L is the standard Laplacian and D is the degree matrix
"""
function find_y_lambda2(A::SparseMatrixCSC)
    L,Dhalf = nLaplacian(A,true)
    n = size(L,1)
    sc = n
    Vl,Vc,convinfo = eigsolve(L + sc*LinearAlgebra.I, 2, :SR; tol = 1e-8, maxiter = 1000, verbosity = 0)
    lam2 = Real(Vl[2])-sc
    x2 = Real.(Vc[2])
    y = Dhalf*x2
    return lam2, y
end

"""
Run a sweepcut procedure on the vector y.
This checks the conductance for each set 

    Si = set of i indices with largest entries in y
    (e.g., if y = [-1 2 9 -2], then S1 = [3]; S2 = [2; 3]; S3 = [1;2;3])

And returns the one with the smallest conductance

Inputs

    A: n x n adjacency matrix for a connected simple graph
    d: degree vector for matrix A
    y: vector of length n on which "sweep" is performed

Outputs

    S:     Si set with minimum conductance
    condS: conducance of S
"""
function sweepcut(A::SparseMatrixCSC,d::Vector,y::Vector)

    p = sortperm(y,rev = true)
    sort_d = d[p]
    bestcond = 1
    bestS = [p[1]]
    volA = sum(d)
    
    S = bestS
    Asrt = A[p,p]               # sorted version of A
    volS = sum(sort_d[1])
    interior_edges = 0
    n = size(A,1)
    for i = 2:n-1
        new_edges = sum(Asrt[1:i-1,i])  # add the new edges
        interior_edges += new_edges     # update interior edge count
        push!(S,p[i])                   # add the new node to S
        volS += sort_d[i]               # increase the denominator
        cutS = volS - 2*interior_edges
        condS = cutS/(min(volS,volA-volS))
        if condS < bestcond
            bestS = copy(S)
            bestcond = condS
        end
    end
    return bestS, bestcond
end

"""
Apply spectral clustering to cluster a simple connected graph.

Input
    
    A: n x n adjacency matrix for a connected simple graph

Output

    lam2:  second smallest eigenvalue of the normalized Laplacian
    y:     vector satisfying lam2 = y'*L*y/(y'*D*y) 
           where L is the standard Laplacian and D is the degree matrix
    S:     Si set with minimum conductance among sweep cuts for y
    condS: conducance of S
"""
function spectral_clustering(A::SparseMatrixCSC)
   
    d = vec(sum(A,dims = 1))
    lam2, y = find_y_lambda2(A)
    S, condS = sweepcut(A,d,y)

    return S, condS, y, lam2
end

"""
Given a graph with adjacency matrix A and a subset of nodes S,
compute the conductance for S

    d = A*ones(n) is the degree matrix of the graph 
    volA = sum(d) is the total volume of the graph
"""
function compute_conductance(A,S,d,volA)
    ES_doubled = sum(A[S,S])
    volS = sum(d[S])
    cutS = volS-ES_doubled
    return cutS/(min(volS,volA-volS))
end