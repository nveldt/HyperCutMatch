## Algorithms for solving the SDP relaxation of hypergraph expansion
using Mosek
using JuMP, MosekTools
using LinearAlgebra

"""
Solve the hypergraph expansion SDP relaxation using Mosek as a backend.
"""
function HypergraphExpansionSDP(H,time_limit=100000,outputflag = true)

    m,n = size(H)
    Elist = incidence2elist(H)

    tic = time()
    # Mosek parameters: https://docs.mosek.com/8.1/capi/param-groups.html#doc-param-groups
    form = 1  # This tells Mosek to focus on the primal, not the dual. There is no automatic dualizer, so choosing form = 2 does nothing, although in theory this means "solve the dual"
    contol = 1e-8 # default dual feasibility tolerance
    model = Model(optimizer_with_attributes(Mosek.Optimizer, "MSK_DPAR_OPTIMIZER_MAX_TIME" => time_limit, "QUIET" => ~outputflag,"INTPNT_CO_TOL_DFEAS" => contol, "MSK_IPAR_INTPNT_SOLVE_FORM" => form))
    @variable(model,X[1:n,1:n],PSD);
    @variable(model, Y[1:m])
    @objective(model, Min, sum(Y[i] for i = 1:m))

    for i = 1:n-2
        for j = i+1:n-1
            for k = j+1:n
                @constraint(model,X[i,i] >= X[i,j] + X[i,k] - X[j,k])
                @constraint(model,X[j,j] >= X[j,i] + X[j,k] - X[i,k])
                @constraint(model,X[k,k] >= X[i,k] + X[j,k] - X[i,j])
            end
        end
    end
    for t = 1:m
        # e = findall(x->x>0,H[t,:])
        e = Elist[t]
        for i = 1:length(e)
            for j = i+1:length(e)
                ii = e[i];
                jj = e[j];
                @constraint(model,Y[t] >= X[ii,ii] + X[jj,jj] - 2*X[ii,jj])
            end
        end
    end

    for i = 1:n
        for j = i+1:n
            @constraint(model,X[i,j] >= 0)
            @constraint(model,X[i,j] <= X[i,i])
            @constraint(model,X[i,j] <= X[j,j])
        end
    end

    @constraint(model, sum(X[i,i] for i = 1:n) == 1)

    for i = 1:n
        @constraint(model,sum(X[:,i]) <= n/2*X[i,i])
    end
    setuptime = time()-tic

    tic = time()
    JuMP.optimize!(model)
    solvertime = time()-tic

    # Status codes: https://jump.dev/MathOptInterface.jl/stable/reference/models/#MathOptInterface.TerminationStatus
    opt_status = termination_status(model) == MOI.OPTIMAL
    time_status = termination_status(model) == MOI.TIME_LIMIT

    if opt_status
        optX = JuMP.value.(X)
        optY = JuMP.value.(Y)
        obj = sum(optY)
        solverstats = [obj, setuptime, solvertime]       
        return opt_status, optX, optY, solverstats 
    else 
        solverstats = [setuptime, solvertime] 
        return opt_status, time_status, 0, solverstats
    end

end

"""
This is an inefficient way to round a set of distance variables X[i,j] 
into a cut with as small of expansion as possible,
but should not be a bottleneck because solving an LP or SDP
to get the X matrix is much more expensive.

For each node, it orders nodes based on how far other nodes are in distance,
and then checks sweep cuts. I.e., for each node u it computes

S = radius(u, r)

for all values of r and then outputs the best set found.

This also works if X[i,j] is a similarity score and X[i,i] is a large nonzero number
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
Given the positive semidefinite matrix X output from the SDP
relaxation, this decomposes X into vectors and computes distances
between points.
"""
function decompose_X(X)
    # decompose X into vectors
    E, V = eigen(X)
    Er = max.(E,0)
    D = sqrt.(Er)
    V = V*Diagonal(D)

    n = size(X,1)
    Dist = zeros(n,n)
    for i = 1:n
        for j = i+1:n
            Dist[i,j] = norm(V[i,:] - V[j,:])
            Dist[j,i] = Dist[i,j]
        end
    end
    return Dist
end