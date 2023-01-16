These are hypergraphs where nodes are amazon products and hyperedges indicate sets of products reviewed by the same user.

The dataset Amazon_9.mat is the combined hypergraph of all products from the 9 smallest product categories. This dataset was previously considered in the paper:

Parameterized Correlation Clustering in Hypergraphs and Bipartite Graphs
Nate Veldt, Anthony Wirth, David F. Gleich
KDD 2020

The file

Data-processing-code/save-amazon-9.jl

is code for loading in this hypergraph and saving the following processed datasets:

1. Amazon9_lcc.mat is the largest connected component of Amazon_9.mat

2. Each individual hypergraphs Am_i_lcc.mat correspond to a subset of products from the same product category. They can be viewed as sub-hypergraps corresponding to the 9 smallest product categories in the larger Amazon hypergraph considered in the following paper:

Minimizing Localized Ratio Cut Objectives in Hypergraphs
Nate Veldt, Austin R. Benson, Jon Kleinberg
KDD 2020

The 'lcc' stands for 'largest connected component'.

Summary statistics
------------------

Here are summary statistics for the amazon-9 datasets:


Max |e|  Avg |e|  n       m      \sum |e|        Hypergraph (product category)
9        7.6     24      396     3025.0          Amazon_Fashion
8        4.1     41      42      174.0   	       Appliances
10       4.2     81      980     4079.0          All_Beauty
31       6.5     148     458     2965.0          Gift_Cards
30       6.6     157     348     2302.0          Magazine_Subscriptions
52       6.5     800     1823    11881.0         Software
90       7.4     1580    3812    28067.0         Luxury_Beauty
183      9.3     4970    14162   131096.0        Prime_Pantry
89       6.5     5333    11025   72009.0         Industrial_and_Scientific




