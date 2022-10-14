using SparseArrays

include("PushRelabelMaxflow.jl")
# include("GurobiMaxflow.jl")
include("pwl_approx.jl")
function TwoNodeEdge(w)
    # a = weight from each node to the source
    # b = weight from one node to the otehr
    #c = weight from each node to the sink
    a = w[1]/2
    b = w[2] - w[1]/2 - w[3]/2
    c = w[3]/2

    return a,b,c
end

function SparseCard(Edges,EdgesW,n,epsilon;refined=true,flowsolver="pushrelabel",gurobioutput=0)

    """
    Wrapper for SparseCard, for solving a Card-DSFM problem.

        refined = true gives better minimal approximate reductions
                (same density of graph, better reductions)
        flowsolver = choose from push-relabel or solving LP with Gurobi
    """

    # Reduce to a graph s-t cut problem
    tic = time()
    A, svec, tvec, worst_eps = CardDSFM_reduction(Edges,EdgesW,n,epsilon,refined)
    reducetime = time()-tic
    N = size(A,1)

    # Solve with a max-flow solver
    if flowsolver == "pushrelabel"

        tic = time()
        F = maxflow(A,svec,tvec,0)
        flowtime = time()-tic
        flowval = F.flowvalue

        # get source set, requires removing source node s = 1
        S_in_A = source_nodes_min(F,1e-10)[2:end].-1

        # ignore auxiliary nodes
        S = intersect(1:n,S_in_A)
        eS = zeros(n)
        eS[S] .= 1
        cutS = eval_fS(Edges,EdgesW,eS)

    else

        Svec, flowval,f,flowtime = Gubori_MaxFlow(A,svec,tvec,gurobioutput)
        eS = Svec[1:n]

    end
    cutval = eval_fS(Edges,EdgesW,eS)

    # The value of the flow in the approximate graph approximates the
    approxcut = flowval
    return eS, cutval, approxcut, flowtime, reducetime, worst_eps
end

function SymmetricCard_reduction(Edges,EdgesW,n,epsilon,returnIJV=true)
    """
    Reduce a hypergraph with symmetric cardinality-based submodular
    splitting functions to a directed graph preserving cuts.

    Edges[i] gives nodes in the ground set of function i in the sum
    EdgesW[i] gives the corresponding CB function penalties
        EdgesW[i][j] gives function evaluation at input (j-1)

    NOTE:

    N = number of objects in the ground set
    epsilon[i] = approximation guarantee for reducing a function/hyperedge of size i
    """
    if length(epsilon) == 1
        K = maximum(length(e) for e in Edges)
        epsilon = epsilon*ones(K)
    end
    N = n
    I = Vector{Int64}()
    J = Vector{Int64}()
    W = Vector{Float64}()
    for i = 1:length(Edges)
        e = Edges[i]
        w = EdgesW[i]
        k = length(e)
        @assert(is_increasing(w))
        @assert(check_cb_submodular(w))
        # @assert(length(w) == floor(Int64,(k)/2)+1)
        # @assert(w[1] == 0)
        if w[1] != 0
            println("For symmetric splitting functions, include the zero: w[1] = w_0 = 0.")
            return
        end
        if length(w) != floor(Int64,(k)/2)+1
            println("Give one penalty for each possible small side of the cut.")
            wrong = length(w)
            right = floor(Int64,(k)/2)+1
            println("length(w) should be $right, but it's $wrong")
            return
        end

        if k == 2
            # just an edge
            push!(I,e[1])
            push!(J,e[2])
            push!(W,w[2])
            push!(I,e[2])
            push!(J,e[1])
            push!(W,w[2])

        elseif k == 3
            # just a triangle
            for ii = 1:3
                for jj = 1:3
                    if ii != jj
                        push!(I,e[ii])
                        push!(J,e[jj])
                        push!(W,w[2]/2)
                    end
                end
            end
        else
            a, b = SymmetricSCB_to_Gadget(w,epsilon[k])

            @assert(minimum(a) > 0)
            @assert(minimum(b) > 0)

            L = length(a)
            for j = 1:L
                # add two new auxiliary nodes for this gadget
                push!(I,N+1)
                push!(J,N+2)
                push!(W,a[j]*b[j])
                for v in e
                    push!(I,v)
                    push!(J,N+1)
                    push!(W,a[j])
                    push!(I,N+2)
                    push!(J,v)
                    push!(W,a[j])
                end
                N += 2
            end
        end
    end

    A = sparse(I,J,W,N,N)
    if returnIJV
        return A,I,J,W
    else
        return A
    end
end

function CardDSFM_reduction(Edges,EdgesW,n,epsilon,refined=true)
    """
    Reduce an instance of Card-DSFM to a graph where the s-t cut
    structure approximately models the cardinality-based decomposable
    submodular function.

    A minimum s-t cut in the graph will then solve the problem.

    Think of the problem as a hypergraph cut problem.

    Edges[i] gives nodes in the ground set of function i in the sum
    EdgesW[i] gives the corresponding CB function penalties
        EdgesW[i][j] gives function evaluation at input (j-1)
    N = number of objects in the ground set
    epsilon[i] = approximation guarantee for reducing a function/hyperedge of size i
    refined = true takes a tiny bit long, but gives better gadget reductions.
        (it finds a better approximation, AND a minimum number of gadgets)
    """
    N = n
    I = Vector{Int64}()
    J = Vector{Int64}()
    W = Vector{Float64}()
    svec = zeros(n)
    tvec = zeros(n)
    if refined
        worst_eps = 0
    else
        worst_eps = maximum(epsilon)
    end
    for i = 1:length(Edges)
        e = Edges[i]
        w = EdgesW[i]
        k = length(e)
        @assert(check_cb_submodular(w))
        if k == 2
            z0,uv,zk = TwoNodeEdge(w)
            svec[e[1]] += z0
            svec[e[2]] += z0
            tvec[e[1]] += zk
            tvec[e[2]] += zk
            push!(I,e[1])
            push!(J,e[2])
            push!(W,uv)
            push!(I,e[2])
            push!(J,e[1])
            push!(W,uv)
        else
            if refined
                # don't just find minimum number of CB-gadgets,
                # find best approx with minimum number of gadgets
                z0, zk, a, b, cgf, eps_val = Refined_SCB_to_CGF(w,epsilon[k])
                if eps_val > worst_eps
                    worst_eps = eps_val
                end
            else
                z0, zk, a, b, cgf = SubCardFun_to_CGF_weights(w,epsilon[k],false)
            end
            @assert(min(z0,zk) > -1e-14)
            z0 = max(z0,0)
            zk = max(zk,0)
            L = length(a)
            if z0 < 0
                @show z0
            end
            if zk < 0
                @show zk
                @show w[end]
            end
            if length(a) > 0 && minimum(a) < 0
                @show minimum(a)
            end
            if length(b) > 0 && minimum(b) < 0
                @show minimum(b)
            end
            for v in e
                svec[v] += z0
                tvec[v] += zk
            end
            for j = 1:L
                for v in e
                    push!(I,v)
                    push!(J,N+j)
                    push!(W,a[j]*(k-b[j]))
                    push!(I,N+j)
                    push!(J,v)
                    push!(W,a[j]*b[j])
                end
            end
            N += L
        end
    end

    A = sparse(I,J,W,N,N)
    svec = [svec;zeros(N-n)]
    tvec = [tvec;zeros(N-n)]
    return A, svec, tvec, worst_eps
end


function eval_fS(Edges::Vector{Vector{Int64}},EdgesW::Vector{Vector{Float64}},eS::Vector)
    """
    Evaluate the cardinality-based decomposable submodular function f at
    set S, where the function is given by ground sets in Edges, each
    with the cardinality-based penalties given in EdgesW.
    """
    cutval = 0
    for i = 1:length(Edges)
        e = Edges[i]
        w = EdgesW[i]
        # t = length(intersect(e,S))
        t = round(Int64,sum(eS[e]))
        cutval += w[t+1]
    end

    return cutval
end


function eval_min_dircut(A,svec,tvec,S,T)
    """
    A is the adjacency matrix between non terminal nodes.
    svec and tvec give terminal edges.
    Fix S to be on the source side, and T to be on the sink side.
    Find the minimum S-T cut with these nodes anchored in this way.
    """
    svec = copy(svec)
    tvec = copy(tvec)
    svec[S] .= 2*sum(A)
    tvec[T] .= 2*sum(A)
    F = maxflow(A,svec,tvec,0)
    S = F.source_nodes[2:end].-1
    cutval = F.cutvalue
    return S, cutval
end
