using MatrixNetworks
using SparseArrays

"""
Take the largest connected component in the star expansion
and use it to get a connected sub-hypergraph.
"""
function largest_cc_star(H)
    m,n = size(H)
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

    return H, d, order
end

"""
Take the largest connected component in the star expansion
and use it to get a connected sub-hypergraph.
"""
function largest_cc_star(H,labels)
    m,n = size(H)
    A = [spzeros(n,n) sparse(H'); H spzeros(m,m)]
    Acc,p = largest_component(A)
    pnodes = p[1:n]
    pedges = p[n+1:n+m]
    H = H[pedges,pnodes]
    Labels = labels[pnodes]
    d = vec(sum(H,dims = 1))
    order = vec(sum(H,dims = 2))
    keep_edges = findall(x->x>1,order)
    H = H[keep_edges,:]
    d = vec(sum(H,dims = 1))
    order = round.(Int64,vec(sum(H,dims = 2)))

    return H, d, order, Labels
end

"""
Iteratively remove a min-degree node and its neighbors until we get to the k-core.

Inefficient code, but not a bottleneck.
"""
function hyperkcore_multi(H,k,verbose = false)
   
    order = vec(round.(Int64,sum(H,dims = 2)))
    keepedges = vec(findall(x->x>1,order))
    H = H[keepedges,:]
    d = vec(round.(Int64,sum(H,dims = 1)))
    m,n = size(H)
    mind = minimum(d)
    nor = n
    it = 1
    while mind < k
        badnodes = vec(findall(x->x==mind,d))
        sumedge = vec(sum(H[:,badnodes],dims = 2))
        badedges = vec(findall(x->x>0,sumedge))
        keepedges = setdiff(collect(1:m),badedges)
        keepnodes = setdiff(collect(1:n),badnodes)
        H = H[keepedges,keepnodes]
        d = vec(round.(Int64,sum(H,dims = 1)))
        m,n = size(H)
        mind = minimum(d)
        if verbose
            println("$it $nor min degree = $mind")
        end
        it += length(badnodes)
    end

    return largest_cc_star(H)

end

"""
This converts a hypergraph in incidence matrix form to hyperedge list form.
Incidence matrix form:  H[e,u] = 1  iff node u is in hyperedge e
Edgelist form: Hyperedges[j] = array of nodes in hyperedge j
nodelist == true: means you actually want to get a map from node IDs to hyperedge ids
"""
function incidence2elist(H::SparseArrays.SparseMatrixCSC{Float64,Int64},nodelist::Bool=false)
    if ~nodelist
        # unless you want the node2edge map, transpose first
        H = SparseArrays.sparse(H')
    end
    rp = H.rowval
    ci = H.colptr
    nz = H.nzval
    Hyperedges = Vector{Vector{Int64}}()
    n,m = size(H)

    for i = 1:m
        startedge = ci[i]
        endedge = ci[i+1]-1
        nodes = rp[startedge:endedge]
        mult = nz[startedge:endedge]
        edge = Vector{Int64}()
        # need to adjust for multiplicities
        for t = 1:length(nodes)
            node = nodes[t]
            for k = 1:mult[t]
                push!(edge,node)
            end
        end

        push!(Hyperedges,edge)
    end
    return Hyperedges
end

"""
Converts a hyperedge list into a binary incidence matrix for the hypergraph.
This is the exact inverse of incidence2elist
"""
function elist2incidence(Hyperedges::Vector{Vector{Int64}}, N::Int64)
    U = Vector{Int64}()
    E = Vector{Int64}()
    M = length(Hyperedges)
    for enum = 1:length(Hyperedges)
        e = Hyperedges[enum]
        for node in e
            push!(U,node)
            push!(E,enum)
        end
    end

    H = SparseArrays.sparse(E,U,ones(length(U)),M,N)
    return H
end


"""
Computes the all-or-nothing hypergraph expansion for S

Takes the hyperedge list and an indicator vector eS as input.
"""
function aon_expansion(Elist,eS)
    n = length(eS)
    S = findall(x->x==eS[1],eS)
    minS = min(length(S), n - length(S))
    cutS = 0
    for t = 1:length(Elist)
        edge = Elist[t]
        s = eS[edge[1]]
        for i = 2:length(edge)
            if eS[edge[i]] != s
                cutS += 1
                break
            end
        end
    end
    sc = cutS*( 1/length(S) + 1/(n-length(S)))
    return cutS/minS, sc
end


"""
This creates a set of symmetric hyperedge splitting penalties, one for each
hyperedge.

EdgesW[j][i] = penalty for splitting hyperedge j such that i-1 nodes are on the small side of the cut

EdgesW[j][1] = 0 always
"""
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

"""
This defines splitting penalties based on the parameterized function
    of Li and Milenkovic (ICML 2018)
"""
function limil_EdgesW(order::Vector{Int64},alpha)
    # fun is a generic function
    EdgesW = Vector{Vector{Float64}}()
    for k in order
        r = floor(Int64,(k)/2)
        w = zeros(r+1)
        thresh = ceil(alpha*k)
        for j = 2:r+1
            # number of nodes on small side of split
            i = j-1 
            w[j] = 1/2 + 1/2*min(1, i/thresh)
        end
        push!(EdgesW,w)
    end
    return EdgesW
end

"""
Delta linear splitting function edges
"""
function Delta_EdgesW(order::Vector{Int64},delta)
    EdgesW = Vector{Vector{Float64}}()
    for k in order
        r = floor(Int64,(k)/2)
        w = zeros(r+1)
        for i = 1:r
            w[i+1] = min(i,delta)
        end
        push!(EdgesW,w)
    end
    return EdgesW
end

"""
alpha-parameterized splitting function from the work of Li and Milenkovic (ICML 2018)
"""
function Alpha_EdgesW(order::Vector{Int64},alpha)
    EdgesW = Vector{Vector{Float64}}()
    for k in order
        r = floor(Int64,(k)/2)
        w = zeros(r+1)
        for i = 1:r
            w[i+1] = 1/2 + 1/2*min(1,i/ceil(alpha*k))
        end
        push!(EdgesW,w)
    end
    return EdgesW
end

"""
For a set S, evaluate generalized hypergraph cut using splitting functions
from EdgesW (symmetric and cardinality-based).
"""
function gen_hypergraph_cut(Edges,EdgesW,S,n)

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


"""
For a set S, evaluate generalized hypergraph ratio cut using splitting functions
from EdgesW (symmetric and cardinality-based), and node-weights from
the vector nodeweights
"""
function gen_ratio_cut(Edges::Vector{Vector{Int64}},EdgesW::Vector{Vector{Float64}},S::Vector{Int64},nodeweights::Vector{Float64},n::Int64,volA::Float64)

    cutS = gen_hypergraph_cut(Edges,EdgesW,S,n)
    volS = sum(nodeweights[S])

    return max(cutS/volS, cutS/(volA-volS))
    
end



"""
Computes the degree of each node in a hypergraph
based on the definition that

d_v = sum_{v: v in e} w_e({v}).

When each hyperedge has a penalty of 1 for splitting off
a single node, then this coincides with the number
of adjacent hyperedges.
"""
function generalized_degree(Edges::Vector{Vector{Int64}},EdgesW::Vector{Vector{Float64}},n::Int64)

    d = zeros(n)

    for t = 1:length(Edges)
        edge = Edges[t]
        weight = EdgesW[t][2]
        for j in edge
            d[j] += weight
        end
    end
    return d
end


## Delta-Linear (tl = thresholded linear) conductance computation.
# e.g. tl_cond(H,S,d,delta,volA,order)
function tl_cond(H::SparseMatrixCSC,S::Vector{Int64},d::Vector{Float64},delta::Float64,volA::Float64,order::Vector{Int64})

    if volA == 0.0
        volA = sum(d)
    end
    n = length(d)
    volS = sum(d[S])
    cut = tl_cut(H,S,delta,order)

    cond = cut/min(volS, volA-volS)

    return cond, volS, cut

end


# Thresholded linear cut value for a set
# calling e.g. tl_cut(H,S,delta,order)
function tl_cut(H::SparseMatrixCSC{Float64,Int64}, S::Vector{Int64}, delta::Float64,order::Vector{Int64})

    # Check the cut
    HS = H[:,S]
    sumHS = sum(HS,dims = 2)  # Count number of S nodes in each hyperedge
    inds = findall(x->x>0,sumHS)    # Extract hyperedges with > 0 nodes from S
    ES = sumHS[inds]
    verts = order[inds]               # Get the size of these hyperedges

    # Find number of nodes on small side of cut
    SmallSide = round.(Int64,min.(ES, verts-ES))
    # Compute the cardinality-based cut score
    cutval = 0.0
    for j = 1:length(SmallSide)
        sm = SmallSide[j]
        if sm > 0
            if sm < delta
                cutval += sm
            else
                cutval += delta
            end
        end
    end

    return cutval
end

