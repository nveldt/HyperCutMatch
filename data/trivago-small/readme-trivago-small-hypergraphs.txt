Each small Trivago hypergraph corresponds to a set of vacation rentals (nodes) that are reviewed by the same user (hyperdges).

Each individual hypergraph is specifically the set of vacation rentals corresponding to the same city, and can be seen as subhypergraphs of the larger TrivagoClickout hypergraph from the paper:

Generative hypergraph clustering: from blockmodels to modularity
Philip S. Chodrow, Nate Veldt, Austin R. Benson
Science Advances, 2021

The original (full hypergraph) dataset is available in text file format in the 'raw-data' folder.

We have *manually* selected a set of cities (node sets) with varying sizes, in order to obtain a set of small benchmark hypergraphs for testing ratio cut objectives.


--------------------
Summary statistics
--------------------

Here is a summary of statistics for the 41 hypergraphs:

Label  max |e|  avg |e|   n      m      sum |e|  Hypergraph (city)
---------------------------------------------------------------------
257      11      3.0     76      43      127.0   Jakarta, Indonesia
1007     12      2.7     92      147     401.0   Angra dos Reis, Brazil
530      9       2.9     95      81      231.0   Villa Gesell, Argentina
50       16      3.1     98      73      225.0   Austin, USA
170      12      3.2     105     139     439.0   Hanover, Germany
515      24      2.9     110     198     579.0   Pigeon Forge, USA
549      10      2.8     116     223     632.0   Recife, Brazil
8        28      3.0     121     101     307.0   Brașov, Romania
2160     11      2.9     125     228     667.0   Puebla, Mexico
2718     16      3.2     130     422     1361.0  Québec City, Canada
742      30      3.3     135     152     506.0   Cairns, Australia
1589     12      3.1     140     234     721.0   Nuremberg, Germany
1019     20      3.2     145     230     740.0   Naha, Japan
993      17      3.1     151     273     852.0   Bruges, Belgium
444      12      2.9     156     286     839.0   Strasbourg, France
1730     17      3.2     160     126     404.0   Hua Hin, Thailand
961      13      3.0     165     473     1399.0  Liverpool, United Kingdom
341      19      3.1     170     295     927.0   Perth, Australia
299      22      3.2     175     203     657.0   San Carlos de Bariloche, Argentina
1365     12      2.7     180     482     1318.0  Durban, South Africa
647      17      3.8     185     282     1058.0  Wrocław, Poland
667      15      3.0     190     233     688.0   Taichung City, Taiwan
651      22      3.0     191     634     1909.0  Sapporo, Japan
765      9       2.6     191     237     611.0   Manila, Philippines
525      22      3.9     195     445     1716.0  Ankara, Turkey
121      13      2.9     200     693     2019.0  Puerto Vallarta, Mexico
1225     16      3.2     206     364     1150.0  Bologna, Italy
442      22      3.1     212     274     860.0   New Orleans, USA
284      16      2.9     216     610     1792.0  Natal, Brazil
1014     10      3.0     221     335     1020.0  Valencia, Spain
155      22      3.4     225     171     584.0   Patong Beach, Thailand
518      19      2.9     229     259     756.0   Cape Town, South Africa
119      16      2.8     235     724     2022.0  Porto de Galinhas, Brazil
88       23      3.1     240     279     863.0   Málaga, Spain
614      24      3.1     245     591     1848.0  Salvador Bahia, Brazil
860      25      3.0     251     635     1874.0  Malacca, Malaysia
1030     21      3.2     257     542     1744.0  Auckland, New Zealand
125      16      3.1     262     910     2833.0  Fukuoka, Japan
559      17      3.0     268     325     987.0   Udaipur, India
379      16      3.1     281     621     1928.0  Cologne, Germany
105      9       3.2     287     244     771.0   Bogotá, Colombia