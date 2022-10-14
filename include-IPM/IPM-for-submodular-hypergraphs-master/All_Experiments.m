addpath('../data/processed')
addpath('Lmats')
addpath('IPM-for-submodular-hypergraphs-master/')
% Load several hypergraphs and convert them to the edge list format
% required by the methods of Li et al.

clear;
%Fix a random seed --- can be changed arbirarily 
random_seed = 1;
rng(random_seed);

names = {'Mushrooms', 'Newsgroups', 'Covertype45'} %, 'Covertype67'};

for i = 1:4
    name = names{i};
    load(strcat(name,'_H.mat'))

    [R,N] = size(H)

    for i = 1:R
        incidence_list{i} = find(H(i,:));
    end

    % Compute the degrees of each vertices
    mu = degree_comp(incidence_list, N, 'full');

    % Control for the outer-loop accuracy: More accurate, Takes more time
    dec_outloop = 1e-1; 
    % Control for the inner-loop accuracy: More accurate, Takes more time
    err_inloop = 1e-2;

    % Run a range of tests for different alphas
    alphas = linspace(0,0.04,11);
    alps = length(alphas);
    IPMs = sparse(N,alps);
    CEs = sparse(N,alps);
    imp_times = zeros(alps,1);
    ce_times = zeros(alps,1);

    for a = 1:11

        alpha = alphas(a);
        load(strcat('Lmats/',name,'_alpha_',num2str(alpha),'.mat'))

        % Set splitting functions
        for i = 1:R
            e = length(incidence_list{i});
            % Compute submodular weights
            parameter_list{i} = card_para([alpha 2], e, 'threshold_bound');
            % Assign types of submodular weights
            submodular_type{i} = 'concave_card';
        end

        [label_inhomo, NCut_inhomo, imp] = run_shp(incidence_list,...
            parameter_list, submodular_type, mu, N, R, ...
            dec_outloop, err_inloop, L, mu);
        IPMs(:,a) = label_inhomo;
        imp_times(a) = imp;

    end
    save(strcat('Output/',name,'_IMP.mat'),'IPM','imp_times')
    
end
    
    