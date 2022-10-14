addpath('../include-IPM/Lmats')
addpath('../include-IPM/IPM-for-submodular-hypergraphs-master/')

clear;

%%
% Fix a random seed --- can be changed arbirarily 

load ../data/trivago-large/trivago_countries_large_summary.mat
random_seed = 1;
rng(random_seed);

inds = [5, 8, 17, 20];
for ii = 1:4

lab = SpecialLabels(inds(ii));    
load(strcat('../data/trivago-large/trivago_countries_',num2str(lab),'_2core.mat'))

[R,N] = size(H);
tic
incidence_list = {}
for i = 1:R
    incidence_list{i} = find(H(i,:));
end
toc
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

    % NEED TO FIRST GENERATE THESE REDUCED GRAPHS BEFORE RUNNING THIS
    % CODE!!
    load(strcat('../data/trivago-large/Lmats/tricount_',num2str(lab),'_',num2str(delta),'_2core_degreewarm.mat'))

    % Set splitting functions
    for i = 1:R
        e = length(incidence_list{i});
        % Compute submodular weights
        parameter_list{i} = delta_linear(delta,e);
        % Assign types of submodular weights
        submodular_type{i} = 'concave_card';
    end

    [eipm, ipmcond, ipmtime] =  run_shp(incidence_list,...
        parameter_list, submodular_type, mu, N, R, ...
        dec_outloop, err_inloop,L,mu);
    sum(eipm)
    save(strcat('../Ex-trivago-large/Output/IPM_tric_',num2str(lab),'_2core_',num2str(delta),'_1_degreewarm.mat'),'ipmcond','ipmtime','eipm')
    
    [ece, ceCond, cetime] = run_ce(incidence_list, parameter_list, submodular_type, N, R, L, mu);
    sum(ece)
    save(strcat('../Ex-trivago-large/Output/CE_tric_',num2str(lab),'_2core_',num2str(delta),'_1_degreewarm.mat'),'ceCond','cetime','ece')
end
end