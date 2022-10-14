function [label_lap, NCut_lap, label_inhomo, NCut_inhomo,ce,imp,L] = demo_fun(...
    incidence_list, parameter_list, submodular_type, mu, N, R, ...
    dec_outloop, err_inloop)

%%% Input 
% incidence_list: incidence matrix, incidence_list{i} contains the list of
% vertices in hyperedge i 

% parameter_list: hyperedge-weight lists, parameter_list{i} contains the list of
% weights of hyperedge i  

% submodular_type: submodular_type{i} indicates What kind of the hyperedge
% weight is using? Currently only support for 'concave_card' 

% mu: degree vector or other positive meansure over vertices
% N: number of vertices
% R: number of hyperedges
% dec_outloop: Control for the outer-loop accuracy: More accurate, Takes more time
% err_inloop: Control for the inner-loop accuracy: More accurate, Takes more time


%%% Output
% label_lap: binary clustering labels obtained by clique projections (Zhou,NIPS2007)
% NCut_lap: NCut obtained by clique projections (Zhou,NIPS2007)
% label_inhomo: binary clustering labels obtained by IPM for submodular hypergraphs (Li,ICML2018)
% NCut_inhomo: NCut obtained by IPM for submodular hypergraphs (Li,ICML2018)
tic
[label_lap, incidence_partition_lap, NCut_lap, v_lap,L] = clique_exp(incidence_list, parameter_list, submodular_type, N, R);
ce = toc

tic
warmstart = v_lap;
[label_inhomo, incidence_partition_inhomo, NCut_inhomo, gap] = submodular_hypergraph_partition(incidence_list, parameter_list, submodular_type, mu, N, R, dec_outloop, err_inloop, warmstart);
imp = toc
