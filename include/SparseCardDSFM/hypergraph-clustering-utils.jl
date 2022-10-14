using HyperModularity

function hypergraph2incidence(H::hypergraph)
    """
    Convert dictionary hypergraph format to binary edge-by-node incidence matrix
    """
    Hyperedges = Vector{Vector{Int64}}()
    weights = Vector{Float64}()
    for key in keys(H.E)
        ed = collect(keys(H.E[key]))
        wt = collect(values(H.E[key]))
        append!(weights,wt)
        append!(Hyperedges,ed)
    end
    N = length(H.D)
    He2n = elist2incidence(Hyperedges,N)
    return He2n, weights
end

function clique_expansion_penalties(k)
    """
    Symmetric clique splitting penalty for hyperedge of size k.
    Only gives first r = floor(k/2) penalties
    """
    @assert(k > 1)
    r = floor(Int64,(k)/2)
    I = collect(0:k)
    w = I .* (k .- I)
    w = w[1:r+1]
    return w
end


function clique_expansion_EdgesW(order::Vector{Int64})
    EdgesW = Vector{Vector{Float64}}()
    for k in order
        weight = 1/(k-1)
        push!(EdgesW,weight*clique_expansion_penalties(k))
    end
    return EdgesW
end

function sqrt_EdgesW(order::Vector{Int64})
    EdgesW = Vector{Vector{Float64}}()
    for k in order
        r = floor(Int64,(k)/2)
        w = sqrt.(collect(0:1:r))
        # w = min.(w,delta)
        push!(EdgesW,w)
    end
    return EdgesW
end

function generic_EdgesW(order::Vector{Int64},fun)
    # fun is a generic function
    EdgesW = Vector{Vector{Float64}}()
    for k in order
        r = floor(Int64,(k)/2)
        w = fun(collect(0:1:r))
        push!(EdgesW,w)
    end
    return EdgesW
end

function clique_expansion_EdgesW(Edges::Vector{Vector{Int64}},weighted=true)
    EdgesW = Vector{Vector{Float64}}()
    for e in Edges
        k = length(e)
        if weighted
            weight = 1/k
        else
            weight = 1
        end
        push!(EdgesW,weight*clique_expansion_penalties(k))
    end
    return EdgesW
end

function delta_linear_EdgesW(Edges::Vector{Vector{Int64}},delta::Float64)
    EdgesW = Vector{Vector{Float64}}()
    for e in Edges
        k = length(e)
        r = floor(Int64,(k)/2)
        w = collect(0:1:r)
        w = min.(w,delta)
        push!(EdgesW,w)
    end
    return EdgesW
end

function gen_hypergraph_cut(Edges,EdgesW,S,n)
    """
    For a set S, evaluate generalized hypergraph cut using splitting functions
    from EdgesW (symmetric and cardinality-based).
    """
    eS = zeros(n)
    eS[S] = ones(length(S))

    # higher-order cut penalties
    cutval=0
    for i = 1:length(Edges)
        e = Edges[i]
        k = length(e)
        w = EdgesW[i]
        t = round(Int64,sum(eS[e]))

        # do this because it's a symmetric splitting function,
        # so w is only half the size
        t = min(t,k-t)
        cutval += w[t+1]
    end

    return cutval

end


# For a set S in a graph with adjacency matrix A, return some information about
# S including its conductance, number of interior edges, volume, and cut.
function set_stats(A::SparseMatrixCSC{Float64,Int64},
    S::Vector{Int64},volA::Float64)

    if volA == 0.0
        volA = sum(A.nzval)
    end

    if length(S) == size(A,1)
        # then we have an indicator vector
        S = findall(x->x!=0,eS)
        AS = A[:,S]
    else
        # then we have a subset
        @assert(minimum(S) >= 1)
        @assert(maximum(S) <= size(A,1))
        AS = A[:,S]
    end

    vol = sum(AS.nzval);
    SAS = AS[S,:]
    edges = sum(SAS.nzval);
    cut = vol-edges

    cond = cut/minimum([vol,volA-vol]);

    return cut, vol, edges, cond

end
