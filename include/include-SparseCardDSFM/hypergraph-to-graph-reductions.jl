include("pwl_approx.jl")
function SymmetricCard_reduction(Edges,EdgesW,n,epsilon,returnIJV=true)
    """
    Reduce a hypergraph with symmetric cardinality-based submodular
    splitting functions to a directed graph preserving cuts.

    Edges[i] gives nodes in the ground set of function i in the sum
    EdgesW[i] gives the corresponding CB function penalties
        EdgesW[i][j] gives function evaluation at input (j-1)

    NOTE:

    N = number of objects in the ground set
    epsilon[i] = approximation guarantee for reducing a function/hyperedge of size i
    """
    if length(epsilon) == 1
        K = maximum(length(e) for e in Edges)
        epsilon = epsilon*ones(K)
    end
    N = n
    I = Vector{Int64}()
    J = Vector{Int64}()
    W = Vector{Float64}()
    for i = 1:length(Edges)
        e = Edges[i]
        w = EdgesW[i]
        k = length(e)
        @assert(is_increasing(w))
        @assert(check_cb_submodular(w))
        # @assert(length(w) == floor(Int64,(k)/2)+1)
        # @assert(w[1] == 0)
        if w[1] != 0
            println("For symmetric splitting functions, include the zero: w[1] = w_0 = 0.")
            return
        end
        if length(w) != floor(Int64,(k)/2)+1
            println("Give one penalty for each possible small side of the cut.")
            wrong = length(w)
            right = floor(Int64,(k)/2)+1
            println("length(w) should be $right, but it's $wrong")
            return
        end

        if k == 2
            # just an edge
            push!(I,e[1])
            push!(J,e[2])
            push!(W,w[2])
            push!(I,e[2])
            push!(J,e[1])
            push!(W,w[2])

        elseif k == 3
            # just a triangle
            for ii = 1:3
                for jj = 1:3
                    if ii != jj
                        push!(I,e[ii])
                        push!(J,e[jj])
                        push!(W,w[2]/2)
                    end
                end
            end
        else
            a, b = SymmetricSCB_to_Gadget(w,epsilon[k])

            @assert(minimum(a) > 0)
            @assert(minimum(b) > 0)

            L = length(a)
            for j = 1:L
                # add two new auxiliary nodes for this gadget
                push!(I,N+1)
                push!(J,N+2)
                push!(W,a[j]*b[j])
                for v in e
                    push!(I,v)
                    push!(J,N+1)
                    push!(W,a[j])
                    push!(I,N+2)
                    push!(J,v)
                    push!(W,a[j])
                end
                N += 2
            end
        end
    end

    A = sparse(I,J,W,N,N)
    if returnIJV
        return A,I,J,W
    else
        return A
    end
end