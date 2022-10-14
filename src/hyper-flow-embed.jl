include("../include/HyperLocal/maxflow.jl")

"""
Given a partition R, node weight vector p, and parameter alpha, find the maximum s-t flow
in the reduced graph graph G(H,R,alpha).

Inputs:
    * A = adjacency matrix for G(H), first n nodes are nodes from V, remaining nodes are auxiliary nodes
    * R = set of nodes in partition
    * Rbar = complement of R in V. I.e., V - R
    * alpha = parameter; searching for set with ratio cut less than alpha
    * nodeweights = node weights vector
    * N = number of nodes in A
    * n = |V|

Output
    * Flow function F
    * Set S that solves the minimum s-t cut problem
"""
function HyperFlowEmbedStep(A::SparseMatrixCSC{Float64,Int64},R::Vector{Int64},Rbar::Vector{Int64},eta::Float64,
    alpha::Float64, N::Int64, nodeweights::Vector{Float64},n::Int64)

        # Directly set up the flow matrix
        sVec = zeros(N)
        tVec = zeros(N)
        sVec[R] .= alpha*nodeweights[R]
        tVec[Rbar] .= alpha*eta*nodeweights[Rbar]
        F = maxflow(A,sVec,tVec,0)
        Src = source_nodes_min(F)[2:end].-1
        S = intersect(1:n,Src)

        return S,F
end

"""
Full HyperFlowEmbed algorithm.

This will also sparsify the graph for you.

    Output:
        * al_best: the best weighted-expansion (e.g., conductance) value returned by the algorithm
        * S_bset: the set that achieves this best ratio cut score
        * F_last: this is the last flow returned, for the iteration in which no improvement was found.
                  This flow can be decomposed to embed a bipartite graph with congestion 1/al into graph A
"""
function HyperFlowEmbed(Edges::Vector{Vector{Int64}},EdgesW::Vector{Vector{Float64}},sparsityeps::Float64,nodeweights::Vector{Float64},R::Vector{Int64})

    # in general, edges weights do not need to be degrees, but we'll call this vector d
    d = nodeweights
    n = length(nodeweights)
    m = length(Edges)
    volA = sum(d)
    volR = sum(d[R])
    @assert(volR <= volA/2)

    eta = volR/(volA - volR)

    A = SymmetricCard_reduction(Edges,EdgesW,n,sparsityeps,false)
    N = round(Int64,size(A,1))

    Rc = setdiff(1:n,R)               # Complement set of R

    cutR = gen_hypergraph_cut(Edges,EdgesW,R,n)
    condR = cutR/(min(volR,volA-volR))

    println("\nRunning HyperFlowOrEmbed")
    println("------------------------------")

    S_best = R
    a_best = condR
    a_old = condR
    still_improving = true
    Iteration = 1
    F_last = 0
    println("Set with pi-expansion $a_best.")
    while still_improving

        still_improving = false

        tic = time()
        S_new, F_last = HyperFlowEmbedStep(A,R,Rc,eta,a_best,N,d,n)
        stime = time()-tic
        cutS = gen_hypergraph_cut(Edges,EdgesW,S_new,n)
        volS = sum(d[S_new])
        condS = cutS/(min(volS,volA-volS))
        a_new = condS

        if a_new < a_old
            still_improving = true
            S_best = S_new
            nS = length(S_best)
            a_old = a_new
            a_best = a_new
            println("Iter $Iteration: |S| = $nS, pi-expansion(S) = $a_new, min-cut took $stime seconds")
        else
            println("Iter $Iteration: Algorithm converged. Last min-cut took $stime sec")
            println("-------------------------------------------------------------------------")
        end
        Iteration += 1
    end

    B = FlowEmbed(F_last,R,a_best,n)
    return S_best, a_best, B

end

"""
Full HyperFlowEmbed algorithm.

This takes a reduced graph as input, so it assumes you've done this step already.

Also, A should not be a sparsified version. This should exactly match Edges and EdgesW.
"""
function HyperFlowEmbed(A::SparseMatrixCSC{Float64,Int64},Edges::Vector{Vector{Int64}},EdgesW::Vector{Vector{Float64}},nodeweights::Vector{Float64},R::Vector{Int64},verbose = false)

    # in general, edges weights do not need to be degrees, but we'll call this vector d
    # d = nodeweights
    n = length(nodeweights)
    volA = sum(nodeweights)
    volR = sum(nodeweights[R])
    @assert(volR <= volA/2)

    eta = volR/(volA - volR)

    N = round(Int64,size(A,1))

    Rc = setdiff(1:n,R)               # Complement set of R

    cutR = gen_hypergraph_cut(Edges,EdgesW,R,n)
    condR = cutR/(min(volR,volA-volR))

    if verbose
        println("\nRunning HyperFlowOrEmbed")
        println("------------------------------")
    end

    S_best = R
    a_best = condR
    a_old = condR
    still_improving = true
    Iteration = 1
    F_last = 0

    while still_improving

        still_improving = false

        tic = time()
        S_new, F_last = HyperFlowEmbedStep(A,R,Rc,eta,a_best,N,nodeweights,n)
        stime = time()-tic
        cutS = gen_hypergraph_cut(Edges,EdgesW,S_new,n)
        volS = sum(nodeweights[S_new])
        condS = cutS/(min(volS,volA-volS))
        a_new = condS

        if a_new < a_old
            still_improving = true
            S_best = S_new
            nS = length(S_best)
            a_old = a_new
            a_best = a_new
            if verbose
                println("Iter $Iteration: |S| = $nS, pi-expansion(S) = $a_new, min-cut took $stime seconds")
            end
        else
            if verbose
                println("Iter $Iteration: Algorithm converged. Last min-cut took $stime sec")
                println("-------------------------------------------------------------------------")

            end
        end
        Iteration += 1
    end

    tic = time()
    B = FlowEmbed(F_last,R,a_best,n)
    toc = time()-tic
    # println("embed time = $toc")

    # Check the flow decomposition
    # tol = 1e-5
    # dBip = vec(sum(B,dims=2))
    # mis = sum(dBip[R]-d[R])
    # mis2 = sum(dBip[Rc]-eta*d[Rc])
    # if abs(mis) >= tol
    #     @show mis
    #     # @show [dBip[R] d[R]]
    #     @show t = maximum(abs.(dBip[R]-d[R]))
    # end
    # @assert(abs(mis) < tol)
    # @assert(abs(mis2) < tol)
    return S_best, a_best, B

end


function gen_hypergraph_cut(Edges::Vector{Vector{Int64}},EdgesW::Vector{Vector{Float64}},S::Vector{Int64},n::Int64)
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

"""
Given a flow that can be embedded in the auxiliary graph for a certain alpha,
    find the bipartite graph that can be embedded with congestion 1/alpha.
"""
function FlowEmbed(Flow::stFlow,R::Vector{Int64},alpha::Float64,n::Int64,tol::Float64 = 1e-10)

    # tol = 1e-6    # Round flow to zero when it is this small

    # The flow matrix Flow.F here is anti-symmetric: F_{ij} = -F_{ij},
    # just because of how the black-box flow algorithm is implemented.
    # We prefer to work with a matrix of nonnegative values where the ij entry is
    # just the nonnegative amount of flow sent from i to j.
    # Fl = Flow.F
    # vv, maxbad = flowconstraints(Flow.F)
    # println("Check worst: $maxbad")
    Fpos_ind = Flow.F .> tol
    If,Jf,Vf = findnz(Fpos_ind)
    F = Flow.F .* Fpos_ind
    dropzeros!(F)

    # Arcs[i] is the list of nodes j such that edge (i,j) has positive flow on it
    Arcs, Fdeg = ConstructAdj(sparse(F'),size(F,1))
    Nf = length(Fdeg)

    # ptr[i] is the pointer to the next node j such that (i,j) has positive flow
    ptr = ones(Int64,Nf)

    # Get the s-t paths
    # Initialize the vectors storing nonzero entries of the bipartite graph we are building.
    BipR = Vector{Int64}()
    BipRbar = Vector{Int64}()
    BipVals = Vector{Float64}()
    # Later we'll do: Bip = sparse(BipR, BipRbar, BipVals, n,n)

    # Initialize path vector and weights
    Path = Vector{Int64}()                  # empty path
    PathWeights = Vector{Float64}()         # weight of each edge in the path
    inpath = zeros(Int64,Nf)                # indicator for whether a node is in the path: 0 means "no", otherwise it is the position of the node in the path
    curr = 1                                # we start from node 1, the source node

    push!(Path,curr)
    inpath[curr] = 1
    pathlength = 1                          # path length in terms of edges
    startcycle = 0

    stillflow = true                        # there is still flow we haven't decomposed
    totalflow = 0                           # total amount of flow we have decomposed already
    error = 0
    # laststep = "nothing"
    while pathlength < Nf && stillflow

                          # let's grow the flow by an edge
        if ptr[curr] > Fdeg[curr]
            # @show ptr[curr], curr, Fdeg[curr]
            # If we got here, then we must have "used" all the flow routes leaving node curr,
            # but it still has flow going into it. There is a numerical error then.
            # @show Path
            # @show PathWeights
            # println("Issue decomposing paths. Try setting tolerance lower")
            
            # Remove the minimum weight edge in this path (will be small), and reset
            minweight, wheremin = findmin(PathWeights)
            println("Removing $minweight flow to maintain flow conservation at nodes")
            for t = 1:length(PathWeights)
                i = Path[t]
                j = Path[t+1]
                inpath[i] = 0
                if F[i,j] - minweight < tol
                    ptr[i] += 1
                    if i == 1 && ptr[1] > Fdeg[1]
                        stillflow = false
                    end
                end
                F[i,j] -= minweight
            end

            # Then we reset so that the path just starts at s again
            pathlength = 1
            Path = [1]
            PathWeights = Vector{Float64}()
            curr = 1
            inpath[curr] = pathlength
            error += minweight
        end
        if stillflow == false
            break
        end

        pathlength += 1 
        next = Arcs[curr][ptr[curr]]        # the next node in the path: follow next neighbor of node curr
        wt = F[curr,next]                   # get the weight of that edge

        # Discard edges that are below tolerance
        # while wt < tol
        #     ptr[curr] += 1
        #     next = Arcs[curr][ptr[curr]]        # the next node in the path: follow next neighbor of node curr
        #     wt = F[curr,next]                   # get the weight of that edge
        # end

        push!(Path,next)                    # grow the path
        push!(PathWeights,wt)               # and set of path weights
        # laststep = "normal"
        if next == Nf                       # we have a path from s=1 to t=Nf

            # from the construction of the auxiliary graph, we know that the second,
            # and second to last node in the path are in R and complement-of-R
            u = Path[2]
            v = Path[end-1]
            # @assert(in(u-1,R))
            # @assert(~in(v-1,R))

            # Find the minimum weight in the path: we'll send this much flow from s to t
            minweight, wheremin = findmin(PathWeights)
            totalflow += minweight

            # This defines an edge in the bipartite graph we are building
            push!(BipR,u-1)
            push!(BipRbar,v-1)
            push!(BipVals,minweight)

            @assert(length(PathWeights) == length(Path)-1)

            # Then substract the minimum flow from every edge in the path
            for t = 1:length(PathWeights)
                i = Path[t]
                j = Path[t+1]
                # ij_wt = PathWeights[t]
                inpath[i] = 0
                # @assert(F[i,j] == ij_wt)
                # @assert(F[i,j] >= minweight)
                if F[i,j] - minweight < tol
                    ptr[i] += 1
                    # if ptr[i] > Fdeg[i]
                    #     @show i, F[i,j] - minweight
                    # end
                    if i == 1 && ptr[1] > Fdeg[1]
                        stillflow = false
                    end
                end
                F[i,j] -= minweight
            end

            # Then we reset so that the path just starts at s again
            next = 1
            pathlength = 1
            Path = [1]
            PathWeights = Vector{Float64}()

            # Before end of loop we'll do:
            #   curr = next;                (so curr = 1)
            #   inpath[curr] = pathlength;  (so inpath[1] = 1)
            #
            # So this starts us our with a fresh path starting from node s=1.
            # laststep = "stpath"
        end

        if inpath[next] > 0                 # then we have found a cycle


            startcycle = inpath[next]       # tells us where in the path it is located
            @assert(startcycle != 1) # we should never have a cycle going back to source node
            # println("Cycle: $(Path[startcycle])")

            # Extract cycle: get subpath that is the cycle
            cyclenodes = Path[startcycle:end]
            cyclewt = PathWeights[startcycle:end]

            # Find minimum weight in the cycle
            minweight, wheremin = findmin(cyclewt)

            # Then substract the minimum flow from every edge in the cycle
            for t = 1:length(cyclenodes)-1
                i = cyclenodes[t]
                j = cyclenodes[t+1]
                # ij_wt = cyclewt[t]
                inpath[i] = 0
                # @assert(F[i,j] == ij_wt)
                # @assert(F[i,j] >= minweight)
                if F[i,j] -  minweight < tol
                    # We used all the flow on the edge (i,j).
                    # Update ptr so that next time we visit node i
                    # we follow an edge with positive flow
                    ptr[i] += 1
                    # if ptr[i] > Fdeg[i]
                    #     println("this one")
                    #     @show i, F[i,j] - minweight
                    # end
                end
                F[i,j] -= minweight
            end

            # Then trim the path back to the part just before we get to the first edge in the cycle.
            # You might think we could save time by only trimming back to the
            # edge that was saturated, but we need to update PathWeights anyways,
            # so this is conceptually easier and has the same complexity.
            Path = copy(Path[1:startcycle])
            PathWeights = copy(PathWeights[1:startcycle-1])
            pathlength = startcycle
            @assert(length(Path) == startcycle)
            # laststep = "cycle"

        end

        curr = next
        inpath[curr] = pathlength
    end

    # Return the bipartite graph we have built
    Bip = sparse(BipR, BipRbar, BipVals, n,n)
    Bip = Bip + Bip'
    # println("$error flow thrown away, $totalflow used")
    return Bip/alpha

end


function flowconstraints(Fl)
    n = size(Fl,1)
    vv = zeros(n-2,1)
    for i = 2:n-1
        vv[i-1] = abs(sum(Fl[:,i]))
    end
    maxbad = maximum(vv)
    return vv, maxbad
end
