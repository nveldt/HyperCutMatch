using Mosek
using JuMP, MosekTools

## Small example

e1 = [1 2 3 4 5];
e2 = [6 7 8 9 10];
e3 = [5 6];
H = zeros(3,10);
H[1,e1] .= 1
H[2,e2] .= 1
H[3,e3] .= 1
H

## Run the sdp relaxation
m,n = size(H)
timelimit = .0001
model = Model(optimizer_with_attributes(Mosek.Optimizer, "QUIET" => false, "INTPNT_CO_TOL_DFEAS" => 1e-7, "MSK_DPAR_OPTIMIZER_MAX_TIME" =>timelimit))
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
    e = findall(x->x>0,H[t,:])
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

JuMP.optimize!(model)
opt_X = JuMP.value.(X)
opt_Y = JuMP.value.(Y)
status = termination_status(model) == MOI.OPTIMAL