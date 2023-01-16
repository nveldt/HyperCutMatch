Demo for IPM methods for submodular hypergraphs
(Reference： Submodular hypergraphs: p-Laplacians, Cheeger inequalities and spectral clustering, )

main.m : main file to run the demo
submodular_hypergraph_partition.cpp : C++ core IPM methods, execute mex 'submodular_hypergraph_partition.cpp' before using 

Currently, the program only supports for the case when the submodular weights only depend
on the cardinality of the input set: e.g., w(S) = w(|S|) = |S||e\S|. So for hyperedge i,
the parameter submodular_type{i} should be assigned as 'concave_card' and parameter_list
should be assigned as parameter_list{i}[k] = w(k) - w(k-1).




