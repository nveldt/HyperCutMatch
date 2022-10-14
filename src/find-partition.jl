## This file has algorithm implementations for finding new partition for the cut matching
# game. There are no hypergraph functions here, because for this step we are building a
# *graph* that can be embedded into the hypergraph.

using ExponentialAction
using StatsBase
using Arpack
using KrylovKit


"""
Simple original FindBisection method of KRV for cut-matching games.

    n = number of nodes in the hypergraph
    Ms = vector of matchings, each one is a sparse adjacency matrix

This is mainly useful for approximating expansion, and not more general ratio cut objectives
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
Uses Arpack for eigenvalue computation.

Given the matrix H output at one interation by the cut-matching procedure
(H is the weighted union of bipartite graphs),
compute the eigenvalue/vector pair for the normalized Laplacian where
we normalized by node weight diagonal matrix Dpi where Dpi[i,i] is the
node weight for node i ("pi" is the node weight function).

If heatkernel = true, this returns the partition using the matrix exponential (i.e., heat kernel approach).

Inputs:

    H = weighted union of bipartite graphs
    Dpi = diagonal matrix of weights to the power -1/2, so Dpi[i,i] = 1/(pi[i]^(-.5)).
    nodeweights = vector of node weights pi
    volA = sum of node weights, volA = sum(nodeweights) = sum(pi)
    heatkernel = true: means use the HK partition

"""
function SpectralOrHeatKernelPartition(H::SparseMatrixCSC{Float64,Int64},Dpi::Diagonal{Float64},nodeweights::Vector{Float64},volA::Float64,heatkernel::Bool=false)
    n = size(H,1)
    
    # Construct the weight-normalized Laplacian for H
    tic = time()
    dH = vec(sum(H,dims = 1))
    DH = Diagonal(dH)
    L = DH-H
    nL = Dpi*L*Dpi
    setup = time() - tic

    # Get eivenvalue for lower bound
    tic = time()
    sc = 100*n
    evals, evecs = eigs(nL + sc*LinearAlgebra.I; nev = 2,which=:SM)
    lam2 = Real(evals[2])-sc
    eigcomp = time()-tic

    # Get partition for next round  
    if heatkernel
        m = randn(n)
        r = m .- mean(m)
        v = exp(-1/2*nL)*r
    else
        v = Real.(evecs[:,2])
    end

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
    getR = time() - tic

    # println("$setup \t $eigcomp \t $getR")

    return R, lam2
end


"""
Uses KrylovKit for eigenvalue computation

Given the matrix H output at one interation by the cut-matching procedure
(H is the weighted union of bipartite graphs),
compute the eigenvalue/vector pair for the normalized Laplacian where
we normalized by node weight diagonal matrix Dpi where Dpi[i,i] is the
node weight for node i ("pi" is the node weight function).

If heatkernel = true, this returns the partition using the matrix exponential (i.e., heat kernel approach).

Inputs:

    H = weighted union of bipartite graphs
    Dpi = diagonal matrix of weights to the power -1/2, so Dpi[i,i] = 1/(pi[i]^(-.5)).
    nodeweights = vector of node weights pi
    volA = sum of node weights, volA = sum(nodeweights) = sum(pi)
    heatkernel = true: means use the HK partition

"""
function SpectralOrHeatKernelPartition_KrylovKit(H::SparseMatrixCSC{Float64,Int64},Dpi::Diagonal{Float64},nodeweights::Vector{Float64},volA::Float64,eigtol::Float64=1e-8,maxits::Int64=1000,heatkernel::Bool=false)
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

    # Get partition for next round  
    if heatkernel
        m = randn(n)
        r = m .- mean(m)
        v = exp(-1/2*nL)*r
    else
        v = Real.(Vc[2])
    end

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

    # println("$setup \t $eigcomp \t $getR")

    return R, lam2
end

"""
Dense eigenvalue computation.
"""
function SpectralOrHeatKernelPartition_dense(H::Matrix{Float64},Dpi::Diagonal{Float64},nodeweights::Vector{Float64},volA::Float64,heatkernel::Bool=false)
    n = size(H,1)
    
    # Construct the weight-normalized Laplacian for H
    dH = vec(sum(H,dims = 1))
    DH = Diagonal(dH)
    L = DH-H
    nL = Matrix(Dpi*L*Dpi)

    # Get eivenvalue for lower bound
    sc = n/2
    nL = Symmetric(nL + sc*LinearAlgebra.I)
    vl, vc = eigen(nL,2:2)
    lam2 = vl[1] - sc

    # Get partition for next round  
    if heatkernel
        m = randn(n)
        r = m .- mean(m)
        v = exp(-1/2*nL)*r
    else
        v = Real.(vec(vc[:,1]))
    end

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
No eigenvalue computation, just heatkernel partitioning.

expaction = true means "use the ExponentialAction function.
"""
function HeatKernelPartition(H::Matrix{Float64},Dpi::Diagonal{Float64},nodeweights::Vector{Float64},volA::Float64,expaction = true)
    n = size(H,1)
    
    # Construct the weight-normalized Laplacian for H
    dH = vec(sum(H,dims = 1))
    DH = Diagonal(dH)
    L = DH-H
    nL = Matrix(Dpi*L*Dpi)

    # Get partition for next round  
    m = randn(n)
    r = m .- mean(m)

    if expaction
        v = expv(-1/2,nL,r)
    else
        # faster package for multiplying by the matrix exponential
        v = exp(-1/2*nL)*r
    end

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

function lowerbound(H,Alphas,nodeweights,dense = true)

    n = size(H,1)
    Dpi = Diagonal(1 ./ sqrt.(nodeweights))

    # Construct the weight-normalized Laplacian for H
    dH = vec(sum(H,dims = 1))
    DH = Diagonal(dH)
    L = DH-H
    if dense
        nL = Matrix(Dpi*L*Dpi)
    else
        nL = Dpi*L*Dpi
    end

    # Get eivenvalue for lower bound
    if dense
        sc = n/2
        nL = Symmetric(nL + sc*LinearAlgebra.I)
        vl, vc = eigen(nL,2:2)
        lam2 = vl[1] - sc
    else
        sc = n
        Vl,Vc,convinfo = eigsolve(nL + sc*LinearAlgebra.I, 2, :SR; tol = 1e-8, maxiter = 1000, verbosity = 0)
        lam2 = Real(Vl[2])-sc
    end

    gamT = sum(1 ./ Alphas)
    almin = minimum(Alphas)
    ApproxFactor = 2*gamT*almin/lam2
    LowerBound = lam2/(2*gamT)

    return ApproxFactor, LowerBound
end