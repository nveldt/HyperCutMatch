function [ece, ceCond, ce] = run_ce(incidence_list, parameter_list, submodular_type, N, R, L,d)
tic
[V,~] = eigs(L,2);
v = d.^(-0.5).* V(:,2)';
eigtime = toc
tic
[ece, ~, ceCond] = optthreshold(incidence_list, parameter_list, submodular_type, d, N, R, v);
ce = toc