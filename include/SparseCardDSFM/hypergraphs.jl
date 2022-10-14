using SparseArrays
using MatrixNetworks

function hypergraph_largest_component(H::SparseMatrixCSC)
    m,n = size(H)
    A = [spzeros(m,m) H; sparse(H') spzeros(n,n)]
    lcc, pcc = largest_component(A)
    p_nodes = pcc[m+1:end]
    p_edges = pcc[1:m]
    H = H[p_edges,p_nodes]
    return H, p_edges, p_nodes
end

# These are some useful functions for converting back and forth between
# different ways to store hypergraphs.
function NeighborList(He2n::SparseArrays.SparseMatrixCSC{Float64,Int64},Hn2e::SparseArrays.SparseMatrixCSC{Float64,Int64})
    """
    NeighborList: given hypergraph (edge x node) indicence and (node x edge)
    incidence, return a list of neighbors for each node.
    """
    m,n = size(He2n)
    rp = He2n.rowval
    ci = He2n.colptr
    Neighbs = Vector{Vector{Int64}}()
    for i = 1:n
        # Get list of edges that node i is in
        Edges = rp[ci[i]:ci[i+1]-1]

        # build a set of neighboring nodes
        NeighbSet = Vector{Int64}()
        for j = 1:length(Edges)
            ej = Edges[j]               # For each node ej adjacent to I...
            nodes = getnodes(Hn2e,ej)   # ...get nodes in edge ej
            append!(NeighbSet,nodes)
        end
        push!(Neighbs, sort(setdiff(unique(NeighbSet),i)))
    end

    return Neighbs
end


function NeighborList(node2edge::Vector{Vector{Int64}},edge2node::Vector{Vector{Int64}})
    """
    NeighborList: given edge to node list and node to edge list, compute
    neighbors
    """
    n = length(node2edge)
    m = length(edge2node)
    Neighbs = Vector{Vector{Int64}}()
    for i = 1:n
        # Get list of edges that node i is in
        Edges = node2edge[i]

        # build a set of neighboring nodes
        NeighbSet = Vector{Int64}()
        for j = 1:length(Edges)
            ej = Edges[j]               # For each node ej adjacent to I...
            nodes = edge2node[ej]       # ...get nodes in edge ej
            append!(NeighbSet,nodes)
        end
        push!(Neighbs, sort(setdiff(unique(NeighbSet),i)))
    end

    return Neighbs
end

function getedges(He2n::SparseArrays.SparseMatrixCSC{Float64,Int64},I::Int64)
    """
    Given the edge-by-node incidence matrix He2n, return the set of hyperedges
    that node I is in.
    """
    first = He2n.colptr[I]
    last = He2n.colptr[I+1]-1
    edges = He2n.rowval[first:last]

    return edges
end

function getnodes(Hn2e::SparseArrays.SparseMatrixCSC{Float64,Int64},J::Int64)
    """
    Given the node-by-edge incidence matrix Hn2e, return the set of nodes
    contained in hyperedge J
    """
    first = Hn2e.colptr[J]
    last = Hn2e.colptr[J+1]-1
    nodes = Hn2e.rowval[first:last]
    mult = Hn2e.nzval[first:last]
    edge = Vector{Int64}()
    # need to adjust for multiplicities
    for t = 1:length(nodes)
        node = nodes[t]
        for k = 1:mult[t]
            push!(edge,node)
        end
    end

    return edge
end

function incidence2elist(H::SparseArrays.SparseMatrixCSC{Float64,Int64},nodelist::Bool=false)
    """
    This converts a hypergraph in incidence matrix form to hyperedge list form.
    Incidence matrix form:  H[e,u] = 1  iff node u is in hyperedge e
    Edgelist form: Hyperedges[j] = array of nodes in hyperedge j

    nodelist == true: means you actually want to get a map from node IDs to hyperedge ids
    """
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

function elist2incidence(Hyperedges::Vector{Vector{Int64}}, N::Int64)
    """
    Converts a hyperedge list into a binary incidence matrix for the hypergraph.

    This is the exact inverse of incidence2elist
    """
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

function RandomSplittingFunctions(Edges,probs)
    ## Generate random splitting functions
    EdgesW = Vector{Vector{Float64}}()
    for edge in Edges
        k = length(edge)
        w = random_scb_function(k,probs[1],probs[2])
        while ~check_cb_submodular(w,0)
            w = random_scb_function(k,probs[1],probs[2])
        end
        push!(EdgesW,w)
    end
    return EdgesW
end


function CliqueExpansion(Edges::Vector{Vector{Int64}},n::Int64,weighted::Bool=true,binary::Bool=false)
    """
    Weighted clique expansion where a hyperedge e is expanded to a
    weighted clique with each edge having weight 1/(|e| - 1)
    """
    I = Int64[]
    J = Int64[]
    V = Float64[]
    for edge in Edges
        k = length(edge)
        for i = 1:k-1
            ei = edge[i]
            for j = i+1:k
                ej = edge[j]
                push!(I,ei)
                push!(J,ej)
                if weighted
                    push!(V, 1/k)
                else
                    push!(V, 1)
                end
            end
        end
    end
    A = SparseArrays.sparse(I,J,V,n,n)
    for i = 1:n
        A[i, i] = 0.0
    end
    SparseArrays.dropzeros!(A)
    A = SparseArrays.sparse(A+A')
    if binary
        I, J, V = SparseArrays.findnz(A)
        A = SparseArrays.sparse(I, J, 1, n, n)
    end
    return A
end
