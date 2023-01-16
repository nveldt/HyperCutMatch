num = 1
load(strcat('../data/amazon-9/Am_',num2str(num),'_lcc.mat'))
fprintf('Amazon dataset number %d \n\n',num)
timelimit = 1800;

[m,n] = size(H);
% First solve the semidefinite program
fprintf('Solving the SDP with cvx \n');
tic
cvx_clear
cvx_solver mosek
% cvx_solver_settings('MSK_DPAR_OPTIMIZER_MAX_TIME', timelimit) 
cvx_begin sdp
    variable X(n,n) symmetric
    variable Y(m)
    minimize (sum(Y))
    subject to

    for i = 1:n-2
        for j = i+1:n-1
            for k = j+1:n
                X(i,i) >= X(i,j) + X(i,k) - X(j,k)
                X(j,j) >= X(j,i) + X(j,k) - X(i,k)
                X(k,k) >= X(i,k) + X(j,k) - X(i,j)
            end
        end
    end
    for t = 1:m
        e = find(H(t,:));
        for i = 1:numel(e)
            for j = i+1:numel(e)
                ii = e(i);
                jj = e(j);
                Y(t) >= X(ii,ii) + X(jj,jj) - 2*X(ii,jj)
            end
        end
    end
    
    for i = 1:n
        for j = i+1:n
            X(i,j) >= 0
            X(i,j) <= X(i,i)
            X(i,j) <= X(j,j)
        end
    end
    
    trace(X) == 1
    
    for i = 1:n
        sum(X(:,i)) <= n/2*X(i,i)
    end
    X == semidefinite(n);
    
cvx_end

timeSDP = toc





