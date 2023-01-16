function [labels, incidence_partition, NCut, v, L] = clique_exp(incidence_list, parameter_list, submodular_type, N, R)
A = zeros(N, N);
degree_vec = degree_comp(incidence_list, N, 'full');
for i = 1:R,
    if (max(incidence_list{i}) > N) | (length(parameter_list{i})==0),
        a = 1;
    end
    A(incidence_list{i}, incidence_list{i}) = A(incidence_list{i}, incidence_list{i}) + parameter_list{i}(1)/length(incidence_list{i});
    i;
end
A = A - diag(diag(A));
d = degree_vec;
D = diag(d.^(-0.5));
L = D*A*D;
[V,sigma] = eigs(L,2);
v = d.^(-0.5).* V(:,2)';
fprintf('About to call CE\n')
[labels, incidence_partition, NCut] = optthreshold(incidence_list, parameter_list, submodular_type, d, N, R, v);

end