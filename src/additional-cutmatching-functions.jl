# This file includes useful functions that were not the most important pieces of code for initial experimental results,
# but are useful alternatives that could be used in the future. 
#
# This includes, for example, functions that use alternative eigenvalue solvers, or avoid eigenvalue computation altogether.

include("HyperCutMatching.jl")

using ExponentialAction
using LinearAlgebra

"""
Simple original FindBisection method of KRV for cut-matching games.

    n = number of nodes in the hypergraph
    Ms = vector of matchings, each one is a sparse adjacency matrix

This is mainly useful for approximating expansion, and not more general ratio cut objectives.

    Pro: simple and likely can be made very fast
    Con: approximation guarantees are not as strong, only designed for simple expansion objective
"""
function FindBisection(Ms::Vector{SparseMatrixCSC{Float64,Int64}},n::Int64)
    t = length(Ms)
    m = randn(n)
    r = m .- mean(m)
    for k = 1:t
        r = 0.5*(Ms[k]*r + r)
    end
    p = sortperm(r)
    half = round(Int64,floor(n/2))
    R = sort(p[1:half])
    return R
end

"""
Cut matching step that uses the heat-kernel partitioning.

Does not include any eigenvalue computation.

This is based loosely on the MATLAB code here:
https://drive.google.com/drive/folders/1RK-Q_8_S6LFxmWC9oRhlZAicNWb-gcax

From the paper:

Practical Nearly-Linear-Time Approximation Algorithms
for Hybrid and Overlapping Graph Clustering
    Ameranis et al., ICML 2022

The cutfind.m function from that code is much more involved. 

WARNING: it is not entirely clear how this application of heat kernel partitioning maps
    to the theory that guarantees an O(log n) approximation algorithm.
    There may also be changes that can be made to make this more practical.
    Nevertheless, this captures the spirit of heat kernel based partitioning
    and even in the worst case is a very effective heuristic for the cut player strategy.
    
    In practice, one can always compute a lower bound certificate on the expander graph that 
    this creates, leading to specific a posteriori approximation guarantees.
"""
function HeatKernelPartition(H::Matrix{Float64},Dpi::Diagonal{Float64},nodeweights::Vector{Float64},volA::Float64)
    n = size(H,1)
    
    # Construct the weight-normalized Laplacian for H
    dH = vec(sum(H,dims = 1))
    DH = Diagonal(dH)
    L = DH-H
    nL = Matrix(Dpi*L*Dpi)

    # Get partition for next round  
    m = randn(n)
    r = m .- mean(m)
    v = expv(-1/2,nL,r)

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

    return R
end


"""
Run the hypergraph cut matching framework using one of three
different eigensolvers (esolver)

    1 = Arpack
    2 = KrylovKitk
    3 = Standard, dense

This will always choose R based on spectral partitioning.
"""
function hypergraphcutmatch(Edges,EdgesW,nodeweights,T,esolver,eigtol = 1e-8, verbose = true)
    
    verystart = time()
    # Initialization
    n = length(nodeweights)
    Dpi = Diagonal(1 ./ sqrt.(nodeweights))
    volA = sum(nodeweights)
    d = nodeweights

    # Initialize containers for output
    LBs = zeros(T-1)
    Lams = zeros(T-1)
    Alphas = zeros(T-1)
    Approx = zeros(T-1)
    Runtimes = zeros(T-1,3)
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
    S,al,B = HyperFlowEmbed(A,Edges,EdgesW,nodeweights,R,false)
    flowtime = time()-tic 
    
    if esolver == 3
        # The third solver uses dense matrices
        H = Matrix(B) 
    else
        # Other two solvers use a sparse matrix
        H = B
    end
    gamT = 1/al
    
    # Next matching player
    tic = time()
    if esolver == 1
        R, lam2 = SpectralOrHeatKernelPartition(H,Dpi,nodeweights,volA)
    elseif esolver == 2
        R, lam2 = SpectralOrHeatKernelPartition_KrylovKit(H,Dpi,nodeweights,volA,1e-4)
    else
        R, lam2 = SpectralOrHeatKernelPartition_dense(H,Dpi,nodeweights,volA)
    end
    eigcomp = time()-tic

    # Store results
    almin = al    
    ApproxFactor = 2*gamT*almin/lam2
    LowerBound = lam2/(2*gamT)
    # LBs[1] = LowerBound
    # Lams[1] = lam2
    # Alphas[1] = al
    # Approx[1] = ApproxFactor
    # Runtimes[1,:] = [flowtime, eigcomp]

    Sbest = S
    
    if verbose 
        println("t \t LowerBound \t ApproxFactor \t alpha \t eig-time \t flow-time")
    end

    for t = 2:T

        # Matching player
        tic = time()
        S,al,B = HyperFlowEmbed(A,Edges,EdgesW,nodeweights,R)
        flowtime = time()-tic

        if al < almin
            almin = al
            Sbest = S
        end
        H += B 
        gamT += 1/al

        # Spectral partition: cut player
        tic = time()
        if esolver == 1
            R, lam2 = SpectralOrHeatKernelPartition(H,Dpi,nodeweights,volA)
        elseif esolver == 2
            R, lam2 = SpectralOrHeatKernelPartition_KrylovKit(H,Dpi,nodeweights,volA,eigtol)
        else
            R, lam2 = SpectralOrHeatKernelPartition_dense(H,Dpi,nodeweights,volA)
        end    
        eigcomp = time()-tic 

        ApproxFactor = 2*gamT*almin/lam2
        LowerBound = lam2/(2*gamT)
        LBs[t-1] = LowerBound
        Lams[t-1] = lam2
        Alphas[t-1] = al
        Approx[t-1] = ApproxFactor
        sofar = time()-verystart
        Runtimes[t-1,:] = [flowtime, eigcomp,sofar]
        if verbose
            println("$t \t $LowerBound \t $ApproxFactor \t $al \t $eigcomp \t $flowtime")
        end
    end
    return LBs, Lams, Runtimes, Alphas, Approx, H, Sbest
end

"""
Run the hypergraph cut matching framework using the heat kernel
method for finding R. No eigenvalues or vectors are computed
"""
function hypergraphcutmatch_heatkernel(Edges,EdgesW,nodeweights,T,verbose = true, expaction = true)
    
    # Initialization
    n = length(nodeweights)
    Dpi = Diagonal(1 ./ sqrt.(nodeweights))
    volA = sum(nodeweights)
    d = nodeweights

    # Initialize containers for output
    Alphas = zeros(T)
    Flowtimes = zeros(T)
    Rtimes = zeros(T)

    # Reduce to a graph that preserves the cut function
    A = SymmetricCard_reduction(Edges,EdgesW,n,0.0,false)

    # Initial cut player
    tic = time()
    R = get_random_initial_R(nodeweights,volA)
    rtime = time()-tic 

    # Matching player
    tic = time()
    S,al,B = HyperFlowEmbed(A,Edges,EdgesW,nodeweights,R,false)
    flowtime = time()-tic 
    H = Matrix(B)
    gamT = 1/al
    

    # Store results
    almin = al    
    t = 1
    Alphas[t] = al
    Rtimes[t] = rtime
    Flowtimes[t] = flowtime

    Sbest = S
    if verbose 
        println("t \t alpha \t cut-play \t match-play")
    end
    for t = 2:T

        # Cut player
        tic = time()
        R = HeatKernelPartition(H,Dpi,d,volA,expaction)
        rtime = time()-tic

        # Matching Player
        tic = time()
        S,al,B = HyperFlowEmbed(A,Edges,EdgesW,nodeweights,R)
        flowtime = time()-tic
        if al < almin
            almin = al
            Sbest = S
        end

        # Embed
        H += B 

        # Store results
        Alphas[t] = al
        Rtimes[t] = rtime
        Flowtimes[t] = flowtime

        if verbose
            println("$t \t $al \t $rtime \t $flowtime")
        end
    end
    @assert(Sbest != 0)
    return Alphas, Rtimes, Flowtimes, H, Sbest
end


"""
Run the hypergraph cut matching framework using the standard findbisection
approach. Designed only for expansion.
"""
function hypergraphcutmatch_findbisection(Edges,EdgesW,T,n,verbose = true)
    
    # Initialization
    nodeweights = ones(n)
    volA = n

    # Initialize containers for output
    Alphas = zeros(T)
    Flowtimes = zeros(T)
    Rtimes = zeros(T)

    # Reduce to a graph that preserves the cut function
    A = SymmetricCard_reduction(Edges,EdgesW,n,0.0,false)

    # Initial cut player
    tic = time()
    R = get_random_initial_R(nodeweights,volA)
    rtime = time()-tic 

    # Matching player
    tic = time()
    S,al,B = HyperFlowEmbed(A,Edges,EdgesW,nodeweights,R,false)
    flowtime = time()-tic  

    # Store results
    almin = al    
    t = 1
    Alphas[t] = al
    Rtimes[t] = rtime
    Flowtimes[t] = flowtime

    Sbest = S
    if verbose 
        println("t \t alpha \t cut-play \t match-play")
    end
    Ms = Vector{SparseMatrixCSC{Float64,Int64}}()
    push!(Ms,B)

    for t = 2:T

        # Cut player
        tic = time()
        R = FindBisection(Ms,n)
        rtime = time()-tic

        # Matching Player
        tic = time()
        S,al,B = HyperFlowEmbed(A,Edges,EdgesW,nodeweights,R)
        flowtime = time()-tic
        if al < almin
            almin = al
            Sbest = S
        end

        # Embed
        push!(Ms,B)

        # Store results
        Alphas[t] = al
        Rtimes[t] = rtime
        Flowtimes[t] = flowtime

        if verbose
            println("$t \t $al \t $rtime \t $flowtime")
        end
    end

    H = Ms[1]
    for j = 2:length(Ms)
        H += Ms[j]
    end

    return Alphas, Rtimes, Flowtimes, H, Sbest
end