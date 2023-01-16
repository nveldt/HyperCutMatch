clear
addpath('../include/IPM-for-submodular-hypergraphs-master/')

%%
%Fix a random seed --- can be changed arbirarily 
random_seed = 1;
rng(random_seed);

gnormstring = '_gnorm_true';
names = {'Newsgroups','Mushrooms','Covertype45', 'Covertype67'};

for i = 3
    name = names{i};
    load(strcat('../data/benchmark-hypergraphs/',name,'_H.mat'))
    
    [R,N] = size(H);

    incidence_list = {};
    for i = 1:R
        incidence_list{i} = find(H(i,:));
    end
    
    % Control for the outer-loop accuracy: More accurate, Takes more time
    dec_outloop = 1e-3; 
    % Control for the inner-loop accuracy: More accurate, Takes more time
    err_inloop = 1e-5;

    % Run a range of tests for different alphas
    alphas = linspace(0,0.04,11);

    for a = 3

        alpha = alphas(a);

        % Set splitting functions
        parameter_list = {};
        submodular_type = {};
        for i = 1:R
            e = length(incidence_list{i});
            % Compute submodular weights
            parameter_list{i} = alpha_linear(alpha, e);
            % Assign types of submodular weights
            submodular_type{i} = 'concave_card';
        end
        
        % Get warmstart and degree weights from CE method
        if a == 1
           load(strcat('Output/CE_',name,'_alpha_0.0',gnormstring,'.mat'))
        else   
            load(strcat('Output/CE_',name,'_alpha_',num2str(alpha),gnormstring,'.mat'))
        end
        warmstart = y;
        mu = nodeweights';
         
        tic
        [eipm, ~, ipmcond, ~] = submodular_hypergraph_partition(incidence_list, parameter_list, submodular_type, mu, N, R, dec_outloop, err_inloop, warmstart);
        ipmtime = toc;
        sum(eipm)        
        save(strcat('Output/IPM_',name,'_alpha_',num2str(alpha),gnormstring,'.mat'),'alphas','ipmcond','eipm','ipmtime')
    end
end

    