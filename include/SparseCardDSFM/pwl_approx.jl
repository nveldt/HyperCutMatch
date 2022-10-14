"""
Functions for finding piecewise linear approximations to concave
    functions at integer points.

    These functions solve the sparsest reduction problem from the paper.
"""

function check_cb_submodular(w::Vector{Float64},tol = 1e-10)
    """
    Given a set of weights defining a cardinality-based function,
    check whether it is submodular.

    w(i) = w_{i-1}
    """
    for i = 2:length(w)-1
        if 2*w[i] - w[i-1] - w[i+1] < -tol
            return false
        end
    end
    return true
end

function check_cb_symmetric(w::Vector{Float64})
    """
    Check whether weights defining a cardinality-based function
    define a symmetric splitting function.

    If there are k weights, this is for a function defined on
    k-1 nodes, since we include w_0.
    """
    k = length(w)
    for i = 1:round(Int64,floor(k/2))
        if w[i] != w[k-i+1]
            return false
        end
    end
    return true
end

function random_CombGad(k::Int64)
    """
    This generates a random combined cardinality-based gadget.

    Returns the parameters as well as the integer
    evaluations.
    """
    J = rand(round(Int64,k/2):k-1)
    a = rand(J)
    b = sort((k-1)*rand(J))
    z0 = 2*rand()
    zk = 2*rand()
    w = CombGad(z0,zk,k,a,b)
    return z0, zk, a, b, w
end

function random_scb_function(k::Int64,prob1=1.0,prob2=1.0)
    """
    Generate a random cardinality-based submodular function.
        k = number of nodes the function is defined on.

    This generates a random combined cardinality-based gadget.

    prob1 = probability that z0 > 0
    prob2 = probability that zk > 0
    """
    J = rand(round(Int64,k/2):k-1)
    a = round.(rand(J),digits = 2)
    while minimum(a) == 0
        a = round.(rand(J),digits = 2)
    end
    b = round.(sort((k-1)*rand(J)),digits = 2)
    while minimum(b) == 0
        b = round.(sort((k-1)*rand(J)),digits = 2)
    end

    if rand() < prob1
        z0 = 2*round.(rand(),digits = 2)
    else
        z0 = 0.0
    end
    if rand() < prob2
        zk = 2*round.(rand(),digits = 2)
    else
        zk = 0.0
    end
    w = CombGad(z0,zk,k,a,b)
    w = w./maximum(w)*3
end

function CombGad(z0,zk,k::Int64,a,b)
    """
    This is a k-node Combined Gadget function.
    Returs the integer values.
    """
    J = length(a)
    @assert(length(a) == length(b))
    @assert(is_increasing(b))

    @assert(b[J] < k)
    @assert(minimum(a) > 0)
    @assert(minimum(b) > 0)
    I = collect(0:k)
    w = z0*(k .- I) + zk*I + sum( a[j]*min.( I*(k-b[j]), (k .- I)*b[j]) for j = 1:J)

    return w
end

function CombGad_func(z0,zk,k::Int64,a,b)
    """
    This is a k-node Combined Gadget function.
    It returns a callable function for any x
    """
    J = length(a)
    if J > 0
        @assert(length(a) == length(b))
        @assert(is_increasing(b))
        @assert(b[J] < k)
        @assert(minimum(a) > 0)
        @assert(minimum(b) > 0)
        f = x->(z0*(k .- x) + zk*x + sum( a[j]*min.( x*(k-b[j]), (k .- x)*b[j]) for j = 1:J))
    else
        f = x->(z0*(k .- x) + zk*x)
    end
    return f
end

function is_increasing(b)
    for i = 1:length(b)-1
        if b[i+1] < b[i]
            return false
        end
    end
    return true
end

is_decreasing(b) = is_increasing(reverse(b))

function plf_to_CombGad(m::Vector{Float64},b::Vector{Float64},k::Int64,f0::Float64)
    """
    Find the Combined Gadget weights for a concave, piecewise linear function f.

    m = slope vector
    b = breakpoints
    k = defines interval [0,k]
    f0 = f(0)

    Output:
    z0, zk, a, b = parameters for Combined Gadget Function
    f = vector of function evaluations at breakpoints.


    More detail
    ------------
    A nonnegative, piecewise linear function f on [0,k] is uniquely determined
    by its breakpoints, slopes, and a single point (0,f(0)).

    Given the slopes m and breakpoints b, and point f(0) of such a function,
    which is additionally concave (meaning m is a decreasing vector).

    We return the parameters z0, zk, a, and b that define the CombGad function
    whose continuous extension is exactly the function f.

    Let J+1 be the number of linear pieces.

    We'll also return the function values at breakpoints: (b_i, f_i)
    """
    @assert(is_decreasing(m))
    J = length(m)-1
    z0 = f0/k
    a = zeros(J)
    f = zeros(J+2)
    f[1] = f0
    bext = [0;b;k]
    for j = 1:J
        a[j] = 1/k*(m[j]-m[j+1])
        f[j+1] = m[j]*(bext[j+1]-bext[j]) + f[j]
    end
    f[J+2] = m[J+1]*(bext[J+2]-bext[J+1])+f[J+1]
    zk = f[end]/k
    return z0, zk, a, b, f[2:end-1]
end


function breakpoints_to_slopes(b,f)
    """
    Given breakpoints and function evaluations at these breakpoints
        for a nonnegative, piecewise linear function on [0,k],
        return the slopes of the curves.

        b[i] = is the ith breakpoint
        f[i] = the piecewise linear function evaluated at b[i]

        Although not technically breakpoints, the first and last entry of b
        should be 0 and k for some integer k, so that we can define the function.
    """
    @assert(length(b) == length(f))
    @assert(minimum(b) == 0)
    @assert(is_increasing(b))
    @assert(b[1] == 0)
    k = b[end]
    @assert(k == round(k))

    J = length(b)-1
    m = zeros(J)
    for i = 1:J
        m[i] = (f[i+1]-f[i])/(b[i+1]-b[i])
    end
    return m

end

# functions for the greedy search
function line_through(x1,y1,x2,y2)
    """
    Return the line through points (x1,y1) and (x2,y2).
    Gives the slope as well as the line as a callable function.
    """
    m = (y2-y1)/(x2-x1)
    inter = y1-m*x1
    l_fun = x->(m*x .+ inter)
    return m, l_fun
end

function line_intersection(m1,m2,x1,x2,y1,y2)
    """
    Given two lines, defined by slopws m1, m2 and
    points (x1,y1), (x2,y2), find the point of intersection.
    """
    if m1 == m2
        println("Lines are parallel")
        return
    end
    b1 = y1 - m1*x1
    b2 = y2 - m2*x2
    x = (b2-b1)/(m1-m2)
    y = m1*x+b1
    return x,y
end


function Refined_SCB_to_CGF(w,epsi=1e-14)
    """
    This allows you to find the CGF with minimum number of linear pieces
    AND the smallest (or approximately the smallest) approximation guarantee
    you can get using that number of pieces.

    It refines the output from one call to the original algoirhtm
    for finding the minimum number of pieces.
    """
    z0, zk, a, b, cgf = SubCardFun_to_CGF_weights(w,epsi,false)
    Jtrue = length(a) + 1
    J = Jtrue
    best_eps = epsi
    # Now find the same number of linear pieces, but better
    # approximations, if possible
    while Jtrue == J && epsi > 1e-10
        epsi = epsi/2
        z0_t, zk_t, a_t, b_t, cgf_t = SubCardFun_to_CGF_weights(w,epsi,false)
        J = length(a_t) + 1
        if J == Jtrue
            z0 = z0_t
            zk = zk_t
            a = a_t
            b = b_t
            cgf = cgf_t
            best_eps = epsi
        end
    end
    return z0, zk, a, b, cgf, best_eps
end

function SubCardFun_to_CGF_weights(w,epsi=1e-14,verbose = false)
    """
    Given function evaluations for a submodular cardinality-based
    function w, return the parameters needed to construct a
    combination of CB-gadgets that approximately model that function.
    """
    # even if the goal is to perfectly fit the function,
    # need to fit it to within a small tolerance for numerical issues
    epsi = max(epsi,1e-14)
    k = length(w)-1
    Lines, Slopes, X, Y, Ranges = get_pwl_approx(w,epsi,verbose)
    z0, zk, a, b, f = Lines_to_CombGadget(Slopes,X,Y,k)
    cgf = CombGad_func(z0,zk,k,a, b)
    return z0, zk, a, b, cgf
end

function Lines_to_CombGadget(Slopes,X,Y,k)
    """
    Given information about a set of lines,
    (namely, the slope and a point for each line)
    return the parameters that describe the piecewise
    linear function as a Combined Gadget Function.
    """
    b1 = Y[1] - Slopes[1]*X[1]
    f0 = b1
    b = zeros(length(X)-1)
    for t = 1:length(X)-1
        x,y = line_intersection(Slopes[t],Slopes[t+1],X[t],X[t+1],Y[t],Y[t+1])
        b[t] = x
    end

    z0, zk, a, b, f = plf_to_CombGad(Slopes,b,k,f0)

    return z0, zk, a, b, f
end

function SymmetricSCB_to_Gadget(w,epsi,verbose = false)
    @assert(is_increasing(w))

    # in the symmetric case we also assume w_0 = 0
    @assert(w[1] == 0)
    r = length(w) - 1
    Lines, M, X, Y, Ranges = get_pwl_approx(w,epsi,false)
    mlast = M[end]

    # In the symmetric case, penalties are increasing
    # and we need to ensure the last linear piece is flat
    if mlast < 0
        # negative sloped lines not needed
        # replace with flat line
        M[end] = 0
        X[end] = r
        Y[end] = w[end]
        # pop!(Lines)
    elseif Ranges[end][1] == r-1 && Ranges[end-1][2] == r-1 && mlast > 0
        # Covered point r-1 twice. Get rid of it
        M[end] = 0
        X[end] = r
        Y[end] = w[end]
        # pop!(Lines)
    elseif mlast > 0
        push!(M,0)
        push!(X,r)
        push!(Y,w[end])
    end
    # push!(Lines, x->w[end]*ones(length(x)))

    b = zeros(length(X)-1)
    for t = 1:length(X)-1
        if M[t] == M[t+1]
            @show r, M
            @show X, Y
            @show mlast
        end
        x,y = line_intersection(M[t],M[t+1],X[t],X[t+1],Y[t],Y[t+1])
        b[t] = x
    end

    a = zeros(length(M)-1)
    for i = 1:length(a)
        a[i] = M[i]- M[i+1]
    end

    return a,b
end

function get_pwl_approx(w,epsi,verbose = false)
    """
    Given penalties w for a submodular cardinality based function,
    find a minimum sized set of linear pieces that provides
    an upper bounding (1+epsi) approximation to w everywhere.
    """
    k = length(w)-1
    @assert(epsi >= 0)
    Lines = Vector()
    Slopes = Vector{Float64}()
    X = Vector{Float64}()
    Y = Vector{Float64}()
    new_i = 1
    Ranges = Vector{Tuple{Int64,Int64}}()

    while new_i <= k+1
        i = new_i
        x2, w2, m, l_fun, new_i = next_line(w,i, epsi; verbose=verbose)

        if i == k || i == k+1
            push!(Ranges,(k-1,k))
        else
            push!(Ranges,(i-1,new_i-2))
        end
        # Save information about the lines
        push!(Lines,l_fun)
        push!(Slopes,m)
        push!(X,x2)
        push!(Y,w2)

    end
    # push!(Ranges,k+1)
    return Lines,Slopes, X, Y, Ranges
end

function next_line(w,i,epsi;verbose = false)
    """
    w is the submodular cardinality-based penalties we're trying to cover.
    i is the next index in w that currently isn't covered by a linear piece.
    Next index i means next integer we need to cover is x[i] = i-1
    """
    k = length(w)-1
    x = collect(0:k)

    # Check for end case 1: x[i] == k-1 or x[i] == k
    if x[i] == k-1 || x[i] == k
        # Add line covering (k-1,w_{k-1}) and (k,w_k)
        xi = x[end-1]
        wi = w[end-1]
        x2 = x[end]
        w2 = w[end]
        m, l_fun = line_through(xi,wi,x2,w2)
        if verbose println("Last line through ($xi,$wi) and ($x2,$w2)") end
        return x2, w2, m, l_fun, k+2
    end

    # Otherwise, start a new line at point (xi,li).
    # This is the worst approx at xi we will accept, maximizing reach.
    xi = x[i]
    li = (1+epsi)*w[i]
    if verbose println("Start a line at ($(xi),$li)") end

    # Visit next point to define a new line
    current = i+1
    x2 = x[current]
    w2 = w[current]
    m, l_fun = line_through(xi,li,x2,w2)
    if verbose println("Candidate line L: through ($x2,$w2)") end

    current += 1

    # Visit points one at a time, greedily finding best cover
    while true
        x_current = x[current]          # integer value we consider
        w_current = w[current]          # true penalty at that value
        l_current = l_fun(x_current)    # approximation provided by current line

        if verbose
            println("Visiting i = $x_current")
            println("(w[$x_current],L[$x_current],(1+eps)w[$x_current]) = ($(w_current),$(l_current),$((1+epsi)*w_current))")
        end

        if l_current < w_current

            if verbose println("Undershot, need new line.") end

            # These are the points defining the new line
            x2 = x_current
            w2 = w_current
            m, l_fun = line_through(xi,li,x2,w2)
            if verbose println("Candidate line L: through ($x2,$w2)") end

            # If this covers last point, we're done
            if x2 == k
                return x2, w2, m, l_fun, k+2
            else
                current += 1
            end

        elseif l_current > (1+epsi)*w_current
            if verbose println("Overshot, return") end

            # we can return the next line, slope m going through point (x2, w2)
            # and also the next point that is now uncovered

            return x2, w2, m, l_fun, current

        elseif l_current <= (1+epsi)*w_current && x_current == k

            # We're done, there's nothing left to cover.
            return x2, w2, m, l_fun, k+2

        else
            # We've covered a point, but there are more points to go.
            current += 1
            x_current = x[current]
        end
    end
end


function clique_round_bound(epsilon)
    """
    Return the upper bound on the number of linear pieces
    needed to bound the clique concave function from [1,k/2]
    for any k.
    """
    return ceil.(log2.(log2.(1 ./ epsilon))) .* ceil.(2*(1 .+ epsilon) ./ sqrt.(2*epsilon))
end
