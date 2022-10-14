include("../include-IPM/ipm-warmstart.jl")

deltas = [1, 1.1, 1.5, 2.0, 5.0, 10.0, 100.0]

M = matread("../data/trivago-large/trivago_countries_large_summary.mat")
Labels = M["SpecialLabels"]
LabelNames = M["LabelNames"]

inds = [5, 8, 17, 20]

for i in inds
    countryind = Labels[i]
    countryname = LabelNames[countryind]
    l = countryind
    mat = matread("../data/trivago-large/trivago_countries_$(l)_2core.mat")
    H = mat["H"]
    m,n = size(H)
    Hyperedges = HyperModularity.incidence2elist(H)
    order = round.(Int64,vec(sum(H,dims = 2)))
    d = vec(sum(H,dims=1))
    println(countryname)
    

    for a = 1:length(deltas)
        delta = deltas[a]
        ## Go through each hyperedge and build a clique expansion matrix
        @time A = CliqueProjectionDelta(Hyperedges,delta,order,n)

        dd = vec(sum(A,dims = 2)) # use degrees from projection, not hypergraph
        D = spdiagm(0=>dd.^(-0.5))
        L = D*A*D
        
        if isinteger(delta)
            del = round(Int64,delta)
            matwrite("../data/trivago-large/Lmats/tricount_$(l)_$(del)_2core_degreewarm.mat", Dict("L"=>L))
        else
            matwrite("../data/trivago-large/Lmats/tricount_$(l)_$(delta)_2core_degreewarm.mat", Dict("L"=>L))
        end
    end

end

