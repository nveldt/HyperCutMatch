# This is the simplest version of the Hypergraph Cut Matching code, boiled down to the most important essentials for experiment reproduction
include("../include/maxflow.jl")
include("../include/include-SparseCardDSFM/hypergraph-to-graph-reductions.jl")
include("hypergraph-helper-functions.jl")

using KrylovKit
using LinearAlgebra

"""
Run the full hypergraph cut matching framework.

This will always choose R based on spectral partitioning, and compute
an approximation bound at each iteration using Theorem 5.3.

Inputs:

    Edges: hyperedge set, Edges[i] is the set of node indices in hyperedge i
    EdgesW: EdgesW[i] is the set of hyperedge splitting penalties for hyperedge W (assumed to be symmetric)
    nodeweights: node weight vector (called pi in the paper, defines pi-expansion)
    T: number of iterations of cut-matching to run
    eigtol: tolerance to use in underlying eigensolver
    verbose: if true, will output updates at each iteration
    returnH: if true, will return the expander graph that is embedded in the hypergraph

Outputs:

    LBs: lower bounds at each iteration
    Lams: second smallest eigenvalue of the Laplacian-like matrix in Theorem 5.3
    Runtimes: keeps track of eigensolver time, flow time, and total-runtime-so-far, for each iteration
    Alphas: best pi-expansion set found in the ith iteration
    Approx: approximation ratio achieved in each iteration 
    Sbest: best pi-expansion set found overall

    optional output:
        H: expander certificate used to get the lower bound (often dense and expensive to store)


"""
function HyperCutMatch(Edges,EdgesW,nodeweights,T,eigtol = 1e-8, verbose = true, returnH = false)
    
    # time keeper
    verystart = time()

    # Initialization
    n = length(nodeweights)
    Dpi = Diagonal(1 ./ sqrt.(nodeweights))
    volA = sum(nodeweights)
    d = nodeweights

    # Initialize containers for output stats at each iteration
    LBs = zeros(T)
    Lams = zeros(T)
    Alphas = zeros(T)
    Approx = zeros(T)
    Runtimes = zeros(T,3)
    volA = sum(nodeweights)

    # Reduce to a graph that preserves the cut function
    tic = time()
    A = SymmetricCard_reduction(Edges,EdgesW,n,1e-10,false)
    redtime = time()- tic
    if verbose
        println("Took $redtime to form reduced graph")
    end

    # Initial cut player
    R = get_random_initial_R(nodeweights,volA)

    # Matching player
    tic = time()
    S,almin,B = HyperFlowEmbed(A,Edges,EdgesW,nodeweights,R,false)
    flowtime = time()-tic 
    
    # update certificate
    H = B
    gamT = 1/almin
    
    # Next matching player
    tic = time()
    R, lam2 = SpectralPartition_CutPlayer(H,Dpi,nodeweights,volA,1e-4)
    eigcomp = time()-tic

    Sbest = S

    # Store results
    ApproxFactor = 2*gamT*almin/lam2
    LowerBound = lam2/(2*gamT)
    Lams[1] = lam2
    Alphas[1] = almin
    Approx[1] = ApproxFactor
    sofar = time()-verystart
    Runtimes[1,:] = [flowtime, eigcomp,sofar]

    # the first time we compute an eigenvalues, 
    # it it small enough that we cannot trust the lower bound 
    # (very sensitive to small errors in the eigenvalue). 
    LBs[1] = 0  

    
    if verbose 
        println("t \t LowerBound \t ApproxFactor \t alpha \t eig-time \t flow-time")
    end

    for t = 2:T

        # Matching player
        tic = time()
        S,al,B = HyperFlowEmbed(A,Edges,EdgesW,nodeweights,R)
        flowtime = time()-tic

        # update the best expansion set found so far
        if al < almin
            almin = al
            Sbest = S
        end
        H += B 
        gamT += 1/al

        # Spectral partition: cut player
        tic = time()
        R, lam2 = SpectralPartition_CutPlayer(H,Dpi,nodeweights,volA,eigtol)   
        eigcomp = time()-tic 

        ApproxFactor = 2*gamT*almin/lam2
        LowerBound = lam2/(2*gamT)
        LBs[t] = LowerBound
        Lams[t] = lam2
        Alphas[t] = al
        Approx[t] = ApproxFactor
        sofar = time()-verystart

        Runtimes[t,:] = [flowtime, eigcomp, sofar]
        if verbose
            println("$t \t $LowerBound \t $ApproxFactor \t $al \t $eigcomp \t $flowtime")
        end
    end

    if returnH
        return LBs, Lams, Runtimes, Alphas, Approx, H, Sbest
    else
        return LBs, Lams, Runtimes, Alphas, Approx, Sbest
    end
end


"""
Cut-player step of the cut matching game. 

Given the matrix H output at one interation by the cut-matching procedure
(H is the weighted union of bipartite graphs),
compute the eigenvalue/vector pair for the normalized Laplacian where
we normalized by node weight diagonal matrix Dpi where Dpi[i,i] is the
node weight for node i ("pi" is the node weight function).

Uses KrylovKit package for eigenvalue computation. Can be made quicker by avoiding 
    an exact eigenvalue computation at each step.

If heatkernel = true, this returns the partition using the matrix exponential (i.e., heat kernel approach).

Inputs:

    H = weighted union of bipartite graphs
    Dpi = diagonal matrix of weights to the power -1/2, so Dpi[i,i] = 1/(pi[i]^(-.5)).
    nodeweights = vector of node weights pi
    volA = sum of node weights, volA = sum(nodeweights) = sum(pi)
    heatkernel = true: means use the HK partition

"""
function SpectralPartition_CutPlayer(H::SparseMatrixCSC{Float64,Int64},Dpi::Diagonal{Float64},nodeweights::Vector{Float64},volA::Float64,eigtol::Float64=1e-8,maxits::Int64=1000)
    n = size(H,1)
    d = nodeweights

    # Construct the weight-normalized Laplacian for H
    dH = vec(sum(H,dims = 1))
    DH = Diagonal(dH)
    L = DH-H
    nL = Dpi*L*Dpi

    # Get eivenvalue for lower bound
    sc = n
    Vl,Vc,convinfo = eigsolve(nL + sc*LinearAlgebra.I, 2, :SR; tol = eigtol, maxiter = maxits, verbosity = 0)
    lam2 = Real(Vl[2])-sc
    v = Real.(Vc[2])

    tic = time()
    p = sortperm(v)
    half = round(Int64,floor(n/2))
    R = p[1:half]
    volR = sum(d[R])

    # Adjust R to make sure it is less than half the graph in volume
    curr = half+1
    while volR < volA/2
        curr_node = p[curr]
        push!(R,curr_node)
        volR += nodeweights[curr_node]
        curr += 1
    end
    while volR > volA/2
        lastr = R[end]
        R = R[1:end-1]
        volR -= nodeweights[lastr]
    end

    return R, lam2
end


"""
This code provides a way to compute the lower bound and a posteriori approximation
guarantee in Theorem 5.3 of the paper.

    H = the expander certificate
    Alphas = the sequence of (node-weighted) expansion values obtained at each iteration
    nodeweights = node weight vector
"""
function lowerbound(H,Alphas,nodeweights)

    n = size(H,1)
    Dpi = Diagonal(1 ./ sqrt.(nodeweights))

    # Construct the weight-normalized Laplacian for H
    dH = vec(sum(H,dims = 1))
    DH = Diagonal(dH)
    L = DH-H
    nL = Dpi*L*Dpi

    # Get eivenvalue for lower bound
    sc = n
    Vl,Vc,convinfo = eigsolve(nL + sc*LinearAlgebra.I, 2, :SR; tol = 1e-8, maxiter = 1000, verbosity = 0)
    lam2 = Real(Vl[2])-sc

    gamT = sum(1 ./ Alphas)
    almin = minimum(Alphas)
    ApproxFactor = 2*gamT*almin/lam2
    LowerBound = lam2/(2*gamT)
    return ApproxFactor, LowerBound
end


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

This will also sparsify the graph for you, if desired

    Input:
        Edges: hyperedge set, Edges[i] is the set of node indices in hyperedge i
        EdgesW: EdgesW[i] is the set of hyperedge splitting penalties for hyperedge W (assumed to be symmetric)
        sparsityeps: parameter controling how good of an approximation you want to the hypergraph cut function.
                     If too small, can sometimes lead to numerical issues. If you want the hypergraph reduction to be exact, should use small error rate like 1e-10 rather than 0
        nodeweights: node weight function (called pi in paper)
        R: set of nodes defining the bisection

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

This A is assumed to perfectly match the splitting penalties in EdgesW (i.e., A should not just be a sparsified approximation for the penalties in EdgesW)
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


"""
Given a flow that can be embedded in the auxiliary graph for a certain alpha,
    find the bipartite graph that can be embedded with congestion 1/alpha.

    Flow = the s-t flow from a previous flow computation
    R = the set of nodes defining the bisection we used when we got Flow
    alpha = the parameter used when we got Flow
    n = number of nodes
    tol = tolerance we use to round flow to zero
    errorcheck: when true, this outputs some additional information on how much flow we have to 
                round away to zero because of small numerical issues
"""
function FlowEmbed(Flow::stFlow,R::Vector{Int64},alpha::Float64,n::Int64,tol::Float64 = 1e-10, errorcheck = false)

    # The flow matrix Flow.F here is anti-symmetric: F_{ij} = -F_{ij},
    # just because of how the black-box flow algorithm is implemented.
    # We prefer to work with a matrix of nonnegative values where the ij entry is
    # just the nonnegative amount of flow sent from i to j.
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
    # Later we'll perform: Bip = sparse(BipR, BipRbar, BipVals, n,n), to get the bipartite graph adjacency marix

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
    error = 0                               # total amount of flow we have to round away to zero because of numerical issues (should be small)
    
    while pathlength < Nf && stillflow

        # In each iteration we grow the flow path by an edge.

        # Confirm that the node "curr" has more outgoing flow
        if ptr[curr] > Fdeg[curr]
            
            # If we are here, then we must have "used" all the flow routes leaving node curr,
            # but it still has flow going into it. There is a numerical rounding issue,
            # some edge weight must be close to zero, so we will round it away, and start a fresh path
            
            # Remove the minimum weight edge in this path (will be small), and reset the search for a path
            minweight, wheremin = findmin(PathWeights)
            if errorcheck
                println("Removing $minweight flow to maintain flow conservation at nodes")
            end

            for t = 1:length(PathWeights)
                i = Path[t]
                j = Path[t+1]
                inpath[i] = 0
                if F[i,j] - minweight < tol
                    ptr[i] += 1
                    if i == 1 && ptr[1] > Fdeg[1]
                        # this means we have fully decompose the flow: there is no more flow leaving s
                        stillflow = false
                    end
                end
                F[i,j] -= minweight
            end
            inpath[Path[end]] = 0 

            # We reset so that the path just starts at source node s again
            pathlength = 1
            Path = [1]
            PathWeights = Vector{Float64}()
            curr = 1
            inpath[curr] = pathlength
            error += minweight
        end

        if stillflow == false
            # we are done if we exhausted all flow out of s
            break
        end

        pathlength += 1 
        next = Arcs[curr][ptr[curr]]        # the next node in the path: follow next neighbor of node curr
        wt = F[curr,next]                   # get the weight of that edge

        push!(Path,next)                    # grow the path
        push!(PathWeights,wt)               # and set of path weights

        if next == Nf                       # we have a path from s=1 to t=Nf

            # from the construction of the auxiliary graph, we know that the second
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
                inpath[i] = 0

                # Some sanity checks:
                # ij_wt = PathWeights[t]
                # @assert(F[i,j] == ij_wt)
                # @assert(F[i,j] >= minweight)

                if F[i,j] - minweight < tol
                    ptr[i] += 1
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

        if inpath[next] > 0            
            # The "next" node is already in the path. Hence we have found a flow cycle

            startcycle = inpath[next]       # tells us where in the path the cycle starts
            @assert(startcycle != 1)        # sanity check: we should never have a cycle going back to source node!

            # println("Cycle: $(Path[startcycle])")

            # Extract cycle: get subpath that is the cycle
            cyclenodes = Path[startcycle:end]
            cyclewt = PathWeights[startcycle:end]

            # if length(cyclewt) == 0
            #     @show Path
            #     @show PathWeights 
            #     ip =  findall(x->x>0,inpath)
            #     @show ip
            #     @show cyclewt
            #     @show startcycle
            #     @show next
            # end

            # Find minimum weight in the cycle
            minweight, wheremin = findmin(cyclewt)

            # Then substract the minimum flow from every edge in the cycle
            for t = 1:length(cyclenodes)-1
                i = cyclenodes[t]
                j = cyclenodes[t+1]
                inpath[i] = 0

                # Some sanity checks:
                # ij_wt = cyclewt[t]
                # @assert(F[i,j] == ij_wt)
                # @assert(F[i,j] >= minweight)

                if F[i,j] -  minweight < tol
                    # We used all the flow on the edge (i,j).
                    # Update ptr so that next time we visit node i
                    # we follow an edge with positive flow
                    ptr[i] += 1
                end
                F[i,j] -= minweight
            end

            # Then trim the path back to the part just before we get to the first edge in the cycle.
            # You might think we could save time by only trimming back to the
            # edge that was saturated, but we need to update PathWeights anyways,
            # so trimming back to the start of the cycle is conceptually easier and has the same time complexity.
            Path = copy(Path[1:startcycle])
            PathWeights = copy(PathWeights[1:startcycle-1])
            pathlength = startcycle
            @assert(length(Path) == startcycle)

        end

        # Get ready for the next loop
        curr = next
        inpath[curr] = pathlength
    end

    # Return the bipartite graph we have built
    Bip = sparse(BipR, BipRbar, BipVals, n,n)
    Bip = Bip + Bip'
    if errorcheck
        println("$error flow thrown away, $totalflow used")
    end
    return Bip/alpha
end

"""
A simple function to check how badly flow constraints are violated.
The flow function works pretty well--flow constraints are typically
satisfied to within a very small tolerance even if the error is not zero.
"""    
function flowconstraints(Fl)
    n = size(Fl,1)
    vv = zeros(n-2,1)
    for i = 2:n-1
        vv[i-1] = abs(sum(Fl[:,i]))
    end
    maxbad = maximum(vv)
    return vv, maxbad
end



"""
Get an intial bisection defined by a set R with 
half the volume of the graph or less.
"""
function get_random_initial_R(nodeweights,volA)
    n = length(nodeweights)
    p = sortperm(randn(n))
    til = round(Int64,n/2)
    R = p[1:til]
    if sum(nodeweights[R]) > volA/2
        R = setdiff(1:n,R)
    end
    return R
end

