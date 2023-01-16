clear
addpath('../include/IPM-for-submodular-hypergraphs-master/')
%%
%Fix a random seed --- can be changed arbirarily 
random_seed = 1;
rng(random_seed);

load ../data/trivago-countries/trivago_countries_large_summary.mat

inds = [5, 8, 17, 20];
gnormstring = '_gnorm_true';

for ii = 1:4
    lab = SpecialLabels(inds(ii));    
    load(strcat('../data/trivago-countries/trivago_countries_',num2str(lab),'_2core.mat'))

    [R,N] = size(H);
    incidence_list = {};
    for i = 1:R
        incidence_list{i} = find(H(i,:));
    end
    % Compute the degrees of each vertices
    mu = degree_comp(incidence_list, N, 'full');

    % Control for the outer-loop accuracy: More accurate, Takes more time
    dec_outloop = 1e-3; 
    % Control for the inner-loop accuracy: More accurate, Takes more time
    err_inloop = 1e-5;

    % Run a range of tests for different alphas
    deltas = [1; 1.1; 1.5; 2; 5; 10; 100];

    for a = 1:numel(deltas)

        delta = deltas(a);

        % Set splitting functions
        for i = 1:R
            e = length(incidence_list{i});
            % Compute submodular weights
            parameter_list{i} = delta_linear(delta,e);
            % Assign types of submodular weights
            submodular_type{i} = 'concave_card';
        end

        if a == 2 || a == 3
            load(strcat('Output/CE_tric_',num2str(lab),'_delta_',num2str(delta),gnormstring,'.mat'))
        else
            load(strcat('Output/CE_tric_',num2str(lab),'_delta_',num2str(delta),'.0',gnormstring,'.mat'))
        end
        warmstart = y;


        tic
        [eipm, ~, ipmcond, ~] = submodular_hypergraph_partition(incidence_list, parameter_list, submodular_type, mu, N, R, dec_outloop, err_inloop, warmstart);
        ipmtime = toc
        
        sum(eipm)
        save(strcat('Output/IPM_tric_',num2str(lab),'_delta_',num2str(delta),gnormstring,'.mat'),'deltas','ipmcond','eipm','ipmtime')
    end
end

    