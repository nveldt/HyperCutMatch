using MatrixNetworks

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


function outputstats(H,S,d,alpha,volA,order,labels)
    n = length(d)
    e1 = ones(n)
    e2 = 2*ones(n)
    e1[S] .= 2
    e2[S] .= 1
    acc1 = LabelAccuracy(e1,labels)
    acc2 = LabelAccuracy(e2,labels)
    acc = max(acc1,acc2)
    err = 1-acc
    pr1, re1, f11 = PRF(T1,S);
    pr2, re2, f12 = PRF(T2,S);
    pr = max(pr1,pr2)
    re = max(re1,re2)
    f1 = max(f11,f12)

    ## Precision, recall, accuracy
    condS, ncutS, volS, cutS = ac_ncut(H,S,d,alpha,volA,order)
    numS = length(S)

    return [numS; volS; condS; ncutS; acc; err; 1-pr; 1-re; f1]
end

"""
Iteratively remove a min-degree node and its neighbors until we get to the k-core.

Inefficient code, not the current bottleneck.
"""
function hyperkcore(H,k,verbose = false)
   
    order = vec(round.(Int64,sum(H,dims = 2)))
    keepedges = vec(findall(x->x>1,order))
    H = H[keepedges,:]
    d = vec(round.(Int64,sum(H,dims = 1)))
    m,n = size(H)
    mind = minimum(d)
    nor = n
    it = 1
    while mind < k
        degv, v = findmin(d)
        # @show degv
        badedges = findall(x->x>0,H[:,v])
        keepedges = setdiff(collect(1:m),badedges)
        keepnodes = setdiff(collect(1:n),v)
        H = H[keepedges,keepnodes]
        d = vec(round.(Int64,sum(H,dims = 1)))
        m,n = size(H)
        mind = minimum(d)
        if verbose
            println("$it $nor min degree = $mind")
        end
        it += 1
    end

    return largest_cc_star(H)

end


"""
Iteratively remove a min-degree node and its neighbors until we get to the k-core.

Inefficient code, not the current bottleneck.
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

## Delta-Linear (thresholded linear) normalized Cut computation.
# e.g. tl_ncut(H,S,d,delta,volA,order)
function tl_ncut(H::SparseMatrixCSC,S::Vector{Int64},d::Vector{Float64},delta::Float64,volA::Float64,order::Vector{Int64})

    if volA == 0.0
        volA = sum(d)
    end
    n = length(d)
    volS = sum(d[S])
    cut = tl_cut(H,S,delta,order)

    cond = cut/min(volS, volA-volS)
    ncut = cut/(volS) + cut/(volA-volS)

    # rncut = round(Int64,ncut)
    # rcut = round(Int64,cut)
    # rcond = round(cond,digits = 4)
    # rvol = round(Int64,volS)

    return cond, ncut, volS, cut

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