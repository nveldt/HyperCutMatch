using StatsBase

include("hyper-flow-embed.jl")
include("find-partition.jl")

include("../include/SparsecardDSFM/hypergraph-clustering-utils.jl")
include("../include/SparseCardDSFM/SparseCard.jl")

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
Run the hypergraph cut matching framework using all three of the eigensolvers,
to compare how and if the results change depending on the situation

    1 = Arpack
    2 = KrylovKitk
    3 = Standard, dense

This will always choose R based on spectral partitioning.
"""
function hypergraphcutmatch_all(Edges,EdgesW,nodeweights,T,verbose = true)
    
    # Initialization
    n = length(nodeweights)
    Dpi = Diagonal(1 ./ sqrt.(nodeweights))
    volA = sum(nodeweights)
    d = nodeweights

    # Initialize containers for output
    LBs = zeros(T)
    Lams = zeros(T,3)
    Alphas = zeros(T)
    Approx = zeros(T)
    Flowtimes = zeros(T)
    Eigtimes = zeros(T,3)
    volA = sum(nodeweights)

    # Reduce to a graph that preserves the cut function
    A = SymmetricCard_reduction(Edges,EdgesW,n,0.0,false)

    # Initial cut player
    tic = time()
    R = get_random_initial_R(nodeweights,volA)
    cuttime = time()-tic

    # Matching player
    tic = time()
    S,al,B = HyperFlowEmbed(A,Edges,EdgesW,nodeweights,R,false)
    flowtime = time()-tic 
    
    H = B
    gamT = 1/al
    almin = al 
    
    # Next cut player
    tic = time()
    R, lam2_1 = SpectralOrHeatKernelPartition(H,Dpi,nodeweights,volA)
    arpacktime = time()-tic

    tic = time()
    R, lam2_2 = SpectralOrHeatKernelPartition_KrylovKit(H,Dpi,nodeweights,volA)
    kk_time = time()-tic

    tic = time()
    R, lam2_3 = SpectralOrHeatKernelPartition_dense(H,Dpi,nodeweights,volA)
    densetime = time()-tic

    # Store results, with lower bounds
    lam2 = lam2_3   
    ApproxFactor = 2*gamT*almin/lam2
    LowerBound = lam2/(2*gamT)
    t = 1
    LBs[t] = LowerBound
    Lams[t,:] = [lam2_1 lam2_2 lam2_3]
    Alphas[t] = al
    Approx[t] = ApproxFactor
    Eigtimes[t,:] = [arpacktime, kk_time, densetime]
    Flowtimes[t] = flowtime

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
        R, lam2_1 = SpectralOrHeatKernelPartition(H,Dpi,nodeweights,volA)
        arpacktime = time()-tic

        tic = time()
        R, lam2_2 = SpectralOrHeatKernelPartition_KrylovKit(H,Dpi,nodeweights,volA)
        kk_time = time()-tic

        tic = time()
        R, lam2_3 = SpectralOrHeatKernelPartition_dense(H,Dpi,nodeweights,volA)
        densetime = time()-tic

        # Store results
        lam2 = lam2_3
        ApproxFactor = 2*gamT*almin/lam2
        LowerBound = lam2/(2*gamT)
        LBs[t] = LowerBound
        Lams[t,:] = [lam2_1 lam2_2 lam2_3]
        Alphas[t] = al
        Approx[t] = ApproxFactor
        Eigtimes[t,:] = [arpacktime, kk_time, densetime]
        Flowtimes[t] = flowtime
        if verbose
            println("$t \t $LowerBound \t $ApproxFactor \t $al \t $eigcomp \t $flowtime")
        end
    end
    @assert(Sbest != 0)
    return LBs, Lams, Eigtimes, Flowtimes, Alphas, Approx, H, Sbest
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