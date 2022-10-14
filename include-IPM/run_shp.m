function [label_inhomo, NCut_inhomo, imp] = run_shp(incidence_list, parameter_list, submodular_type, mu, N, R, ...
    dec_outloop, err_inloop, L,d)

tic
% We'll load in this L from outside, rather than computing it fresh
% A = sparse(N, N);
% degree_vec = degree_comp(incidence_list, N, 'full');
% for i = 1:R
%     A(incidence_list{i}, incidence_list{i}) = A(incidence_list{i}, incidence_list{i}) + parameter_list{i}(1)/length(incidence_list{i});
% end
% fprintf('Done forming A\n',ws)
% A = A - diag(diag(A));
% d = degree_vec;
% D = diag(d.^(-0.5));
% L = D*A*D;
[V,~] = eigs(L,2);
v = d.^(-0.5).* V(:,2)';
warmstart = v;
ws = toc;
fprintf('%f seconds to warm start\n',ws)
tic
[label_inhomo, incidence_partition_inhomo, NCut_inhomo, gap] = submodular_hypergraph_partition(incidence_list, parameter_list, submodular_type, mu, N, R, dec_outloop, err_inloop, warmstart);
imp = toc