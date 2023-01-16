clear;
%Fix a random seed --- can be changed arbirarily 
random_seed = 1;
rng(random_seed);
% generate data
% Number of vertices
N = 100;
% Number of hyperedges
R = 100;
% The size of hyperedge --- Demo only for uniform hypergraphs
edgesize = 10;
for i = 1:R
    % Incidence matrix
    incidence_list{i} = randperm(N, edgesize);
end
% Compute the degrees of each vertices
mu = degree_comp(incidence_list, N, 'full');

% Parameter to control submodular weights over hyperedges --- Currently
% only support the formula in the Section Experiment of the paper
% Submodular hypergraphs ... [ICML2018]
alpha = 0.5;
for i = 1:R
    e = length(incidence_list{i});
    % Compute submodular weights
    parameter_list{i} = card_para([alpha 2], e, 'threshold_bound');
    % Assign types of submodular weights
    submodular_type{i} = 'concave_card';
end

%%
% Control for the outer-loop accuracy: More accurate, Takes more time
dec_outloop = 1e-3; 
% Control for the inner-loop accuracy: More accurate, Takes more time
err_inloop = 1e-6;

[label_lap, NCut_lap, label_inhomo, NCut_inhomo] = demo_fun(...
    incidence_list, parameter_list, submodular_type, mu, N, R,...
    dec_outloop, err_inloop);

fprintf(' random seed:%d \n NCut by Clique Expansion:%f \n NCut by IMP for submodular hypergraphs:%f\n', random_seed, NCut_lap, NCut_inhomo);



