## Implement a Gurobi solver for finding a lower bound on hypergraph all-or-nothing sparsest cut

include("../include/SparsecardDSFM/hypergraph-clustering-utils.jl")
using Gurobi
using JuMP
using LinearAlgebra
gurobi_env = Gurobi.Env()

# Use Gurobi with JuMP interface to get a lower bound for sparsest cut in hypergraphs
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
        return false, termination_status(ml), 0, solverstats
    end

end


"""
This is an inefficient way to round the LP relaxation
into a cut with as small of expansion as possible,
but should not be a bottleneck because solving the LP
is much more expensive.
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
With an indicator vector input eS and hyperedge list Elist,
compute the all-or-nothing hypergraph expansion for S
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


function find_violations!(D::Matrix{Float64}, violations::Vector{Tuple{Int,Int,Int}})
    n = size(D,1)
  
    # We only need this satisfied to within a given tolerance, since the
    # optimization software will only solve it to within a certain tolerance
    # anyways. This can be tweaked if necessary.
    epsi = 1e-8
    @inbounds for i = 1:n-2
         for j = i+1:n-1
            a = D[j,i]
             for k = j+1:n
                b = D[k,i]
                c = D[k,j]
          if a - b > epsi && a - c > epsi && a-b-c > epsi
              push!(violations, (i,j,k))
                  # @constraint(m, x[i,j] - x[i,k] - x[j,k] <= 0)
          end
          if b - a > epsi && b - c > epsi && b-a-c > epsi
              push!(violations, (i,k,j))
              # @constraint(m, x[i,k] - x[i,j] - x[j,k] <= 0)
          end
  
          if c - a > epsi && c-b>epsi && c-a-b > epsi
              push!(violations, (j,k,i))
              # @constraint(m, x[j,k] - x[i,k] - x[i,j] <= 0)
          end
        end
      end
    end
end


# Subroutine in the rounding scheme
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

# Apply the randomized rounding procedure many times
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