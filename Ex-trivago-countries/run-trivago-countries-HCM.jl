using MAT
using SparseArrays
using MatrixNetworks
using LinearAlgebra

## Code
include("../src/HyperCutMatching.jl");

# Load dataset
M = matread("../data/trivago-countries/trivago_countries_large_summary.mat")
Labels = M["SpecialLabels"]
LabelNames = M["LabelNames"]

inds = [5, 8, 17, 20]

for jj = 1:5
    for i in inds
        countryind = Labels[i]
        countryname = LabelNames[countryind]
        l = countryind
        mat = matread("../data/trivago-countries/trivago_countries_$(l)_2core.mat")
        H = mat["H"] 
        m,n = size(H)
        order =  vec(sum(H,dims = 2))
        order = round.(Int64,order)
        d = vec(sum(H,dims = 1))
        @assert(length(d) == n)
        @assert(length(order) == m)
        mu = sum(order)
        maxr = round(Int64,maximum(order))
        meanr = round(sum(order)/length(order),digits = 1)
        println("$countryname: \t  $n \t  $m \t  $maxr \t $meanr \t $mu")
        volA = sum(d)
        Edges = incidence2elist(H);

        # Run Hypergraph Cut Matching
        deltas = [1 1.1 1.5 2 5 10 100]
        verbose = false
        returnH = true
        eigtol = 1e-4           # tolerance for underlying eigenvalue solver

        for a = 1:length(deltas)
            delta = deltas[a]
            EdgesW = Delta_EdgesW(order,delta)     # hyperedge splitting functions
            nodeweights = d                        # node weights (degrees, for conductance)
            T = 10*round(Int64,log2(n))            # number of cut matching iterations
            tic = time()
            LBs, Lams, RuntimesHCM, Conds, Approx, H, S =  HyperCutMatch(Edges,EdgesW,nodeweights,T,eigtol,verbose,returnH);
            solvetime = time()-tic
            ApproxFac, Lbound = lowerbound(H,Conds,nodeweights)

            ratio = round(ApproxFac,digits = 3)
            lb = round(Lbound,digits=3)
            bestcond = round(minimum(Conds),digits = 3)
            numS = min(length(S), n-length(S))
            matwrite("Output/HCM_tric_$(l)_delta_$(delta)_$(jj).mat",Dict("ApproxFac"=>ApproxFac,"Lbound"=>Lbound,"LBs"=>LBs, "Lams"=>Lams, "RuntimesHCM"=>RuntimesHCM,"Conds"=>Conds,"Approx"=>Approx,"S"=>S,"solvetime"=>solvetime))
            println("\t$delta: $lb $bestcond $ratio $solvetime $numS $jj")
        end
    end
end

