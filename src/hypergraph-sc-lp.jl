## Implement a Gurobi solver for finding a lower bound on hypergraph all-or-nothing sparsest cut
# ENV["GUROBI_HOME"] = "/Library/gurobi912/mac64/"
using Gurobi
using JuMP
using LinearAlgebra
gurobi_env = Gurobi.Env()

include("hypergraph-helper-functions.jl")

"""
Uses Gurobi with JuMP interface to get a lower bound for sparsest cut in hypergraphs
"""
function HypergraphLeightonRao(H,time_limit=100000,outputflag = true)

    m,n = size(H)
    Elist = incidence2elist(H)
    m = length(Elist)

    # Create a model that will use Gurobi to solve
    # m = Model(solver=GurobiSolver(OutputFlag = 1,TimeLimit = time_limit, FeasibilityTol = FeasTol, Method = SolverMethod, Crossover = CrossoverStrategy, LogFile = OutputFile))
    tic = time()
    ml = Model(optimizer_with_attributes(() -> Gurobi.Optimizer(gurobi_env), "OutputFlag" => outputflag, "TimeLimit" => time_limit))

    @variable(ml, x[1:n,1:n])
    @variable(ml, y[1:m])

    # Minimize the sum of hyperedge variables
    @objective(ml, Min, sum(y[i] for i = 1:m))

    # Constraint 0 <= x_ij for all node pairs (i,j)
    for i = 1:n-1
        for j = i+1:n
            @constraint(ml,x[i,j] >= 0)
        end
    end

    # Triangle inequality contraint
    for i = 1:n-2
        for j = i+1:n-1
            for k = j+1:n
                @constraint(ml,x[min(j,k),max(j,k)] - x[min(i,k),max(i,k)] - x[min(i,j),max(i,j)] <= 0)
                @constraint(ml,-x[min(j,k),max(j,k)] + x[min(i,k),max(i,k)] - x[min(i,j),max(i,j)] <= 0)
                @constraint(ml,-x[min(j,k),max(j,k)] - x[min(i,k),max(i,k)] + x[min(i,j),max(i,j)] <= 0)
            end
        end
    end
    
    # sum constraint
    @constraint(ml, sum(x[i,j] for i = 1:n-1 for j = i+1:n) == n)

    for t = 1:m
        edge = Elist[t]
        for ii = 1:length(edge)-1
            for jj = ii+1:length(edge)
                ei = edge[ii]
                ej = edge[jj]
                @constraint(ml,y[t] >= x[min(ei,ej),max(ei,ej)])
            end
        end
    end
    setuptime = time()-tic

    tic = time()
    JuMP.optimize!(ml)
    solvertime = time()-tic

    if has_values(ml)
        X = JuMP.value.(x)
        Y = JuMP.value.(y)
        obj = sum(Y)
        solverstats = [obj, setuptime, solvertime]       
        return true, X, Y, solverstats 
    else 
        solverstats = [setuptime, solvertime] 
        @show termination_status(ml)
        return false, termination_status(ml), 0, solverstats
    end

end


"""
This is an inefficient way to round the LP relaxation
into a cut with as small of expansion as possible,
but should not be a bottleneck because solving the LP
is much more expensive.

For each node, it orders nodes based on how far other nodes are in distance,
and then checks sweep cuts. I.e., for each node u it computes

S = radius(u, r)

for all values of r and then outputs the best set found.
"""
function round_hyperLR(X,Elist)
    n = size(X,1)
    Is, Js, Vs = findnz(sparse(X))
    if maximum(Vs) - minimum(Vs) < 1e-12
        x = sparse(Is,Js,1.0,n,n)
        c = extractClustering(x)
        expS, scS = aon_expansion(Elist,c)
        return c, expS
    end

    bpoints = unique(round.(Vs,digits = 4))
    Rx = round.(X + X',digits = 4)
    bestExp = Inf
    bestS = 0

    # do a threshold cut (i.e., sweep cut)
    # for all starting nodes
    for i = 1:n
        dvec = Rx[:,i]
        for b in bpoints
            S = findall(x->x<b,dvec)
            c = zeros(n)
            c[S] .= 1
            expS, scS = aon_expansion(Elist,c)
            if expS < bestExp
                bestExp = expS
                bestS = c
            end
        end
    end

    return bestS, bestExp
end

"""
Avoids setting up the problem using JuMP. JuMP can make things slower.
"""
function HypergraphLR_noJuMP(H,timelimit=100000,outputflag = true)
    m,n = size(H)
    Elist = incidence2elist(H)
    M = length(Elist)
    
    # build map of variables and the objective
    
    # The first m variables correspond to the first m edges
    
    tic = time()
    vmap = -ones(Int,n,n)
    obj = Float64[]
    for t = 1:M
        push!(obj, 1)
    end

    # map from node pairs to their variable index
    nvars = M
    for j=1:n
        for i=1:j-1
            nvars += 1
            vmap[i,j] = nvars-1
            vmap[j,i] = nvars-1
            push!(obj, 0)
        end
    end
    vtypes = repeat(GRB_CONTINUOUS, nvars)
    aptr = Ref{Ptr{Cvoid}}()
    err = GRBnewmodel(gurobi_env, aptr, "HyperLR", nvars, obj, C_NULL, C_NULL, vtypes, C_NULL)
    m = aptr[]
    GRBsetdblparam(GRBgetenv(m), GRB_DBL_PAR_TIMELIMIT, timelimit)
    GRBsetintparam(GRBgetenv(m), "OutputFlag",outputflag)

    # First set of contraints: x[e] >= x[u,v] for each u,v \in e
    for t = 1:M
        edge = sort(Elist[t])
        cind = Int32[0,0]
        cval = Float64[0,0]

        cind[1] = t-1  # indices start from 0, not 1
        cval[1] = -1
        cval[2] = 1
        for ii = 1:length(edge)
            for jj = ii+1:length(edge)
                i = edge[ii]
                j = edge[jj]
                cind[2] = vmap[i,j]

                # x[uv] - x[e] <= 0 for u,v \in e
                error = GRBaddconstr(m, 2, cind, cval, GRB_LESS_EQUAL, 0.0, C_NULL)
            end
        end
    end

    # Second set of constraints: triangle inequalities
    cind = Int32[0,0,0]
    cval = Float64[0,0,0]
    for i = 1:n-2
        for j = i+1:n-1
            for k = j+1:n
                cind[1] = vmap[j,k]
                cind[2] = vmap[i,k]
                cind[3] = vmap[i,j]
                cval[1] = 1
                cval[2] = -1
                cval[3] = -1
                error = GRBaddconstr(m, 3, cind, cval, GRB_LESS_EQUAL, 0.0, C_NULL)

                cval[1] =-1
                cval[2] = 1
                cval[3] = -1
                error = GRBaddconstr(m, 3, cind, cval, GRB_LESS_EQUAL, 0.0, C_NULL)

                cval[1] = -1
                cval[2] = -1
                cval[3] = 1
                error = GRBaddconstr(m, 3, cind, cval, GRB_LESS_EQUAL, 0.0, C_NULL)

            end
        end
    end

    # Final constraint: sum x[ij] = n
    N = Int64(n*(n-1)/2)
    cval = ones(N)
    cind = Int32.(collect(M:M+N-1))
    error = GRBaddconstr(m, N, cind, cval, GRB_EQUAL, n, C_NULL)
    setuptime = time()-tic

    tic = time()
    GRBoptimize(m)
    solvertime = time()-tic

    stat = Ref{Int32}(0)
    GRBgetintattr(m, GRB_INT_ATTR_STATUS, stat)
    # Status codes: https://www.gurobi.com/documentation/9.5/refman/optimization_status_codes.html#sec:StatusCodes
    status = stat[]
    # println("Status = $status")
    if status == 2
        optimal = true
    else
        optimal = false
    end

    if optimal

        robj = Ref{Float64}(0.0)
        GRBgetdblattr(m, GRB_DBL_ATTR_OBJVAL, robj)
        obj = robj[]
        soln = zeros(nvars)
        GRBgetdblattrarray(m, GRB_DBL_ATTR_X, 0, nvars, soln)
        X = zeros(n,n)
        for j=1:n
            for i=1:j-1
                X[i,j] = soln[vmap[i,j]+1]
                X[j,i] = soln[vmap[i,j]+1]
            end
        end
        Y = zeros(M)
        for t = 1:M
            Y[t] = soln[t]
        end
        obj = sum(Y)
        solverstats = [obj, setuptime, solvertime] 
        GRBfreemodel(m)       
        return true, X, Y, solverstats 
    else 
        solverstats = [setuptime, solvertime] 
        GRBfreemodel(m) 
        return false, status, 0, solverstats
    end

end

"""
If x is a 0,1 matrix satisfying the triangle inequality, then we
can extract a clustering from it.
"""
function extractClustering(x)
    # Assuming the triangle inequality results work, we just have to go through
    # each row of x, and if an element hasn't been placed in a cluster yet, it
    # it starts its own new one and adds everything that is with it.
    n = size(x,1)
    NotClustered = fill(true,n)
    c = zeros(n)
    clusnum = 0
    for i = 1:n
        if NotClustered[i]
            for j = i:n
                if x[i,j] < .01 # must be zero, they are clustered together
                    c[j] = clusnum
                    NotClustered[j] = false;
                end
            end
            clusnum +=1
        end
    end
    return round.(Int64,c)
end


# Subroutine in the rounding scheme:
# See https://lucatrevisan.github.io/expanders2016/lecture10.pdf
# for details
function onedim_embedding(X)
    @assert(issymmetric(X))
    n = size(X,1)
    K = round(Int64,log2(n))
    t = rand(1:K)
    p = 1/2^t
    eA = sprand(Bool,n,p)
    while sum(eA) == 0
        eA = sprand(Bool,n,p)
    end
    A,V = findnz(eA)
    DA = X[A,:]
    fA = zeros(n)
    for i = 1:n
        fA[i] = minimum(DA[:,i])
    end
    return fA
end

# Sweep cut procedure (inefficient, but not the bottleneck)
function sweepcut(f,Elist)
    n = length(f)
    bestExp = Inf
    bestS = 0
    fp = round.(f,digits = 4)
    bpoints = sort(unique(fp))

    # threshold/sweep cut
    for b in bpoints
        S = findall(x->x<b,fp)
        c = zeros(n)
        c[S] .= 1
        expS, scS = aon_expansion(Elist,c)
        if expS < bestExp
            bestExp = expS
            bestS = c
        end
    end


    return bestS, bestExp
end

"""
Apply the randomized rounding procedure for the Leighton-Rao-based LP relaxation, many times.
"""
function round_LP_many(X,m,Elist)
    n = size(X,1)
    bestExp = Inf
    bestS = 0

    for b = 1:m
        fA = onedim_embedding(X)
        c, expS = sweepcut(fA,Elist)
        if expS < bestExp
            bestExp = expS
            bestS = c
        end
    end
    return bestS, bestExp
end