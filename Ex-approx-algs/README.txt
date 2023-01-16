Code for running approximation algorithms on benchmark graphs

------------------------------------------------
SDP Relaxation for O(sqrt(log n)) approximation
------------------------------------------------

run_sdp_amazon.m
run_sdp_trivago.m
run_sdp_mathoverflow.m

These are all for running the SDP relaxation in CVX in Matlab, as given by the code in hypergraphexpSDP.m

You can also solve the SDP with Julia as a front end (calling Mosek solver underneath). Unlike CVX, this does not automatically dualize the problem first, so it crashes even on smaller hypergraphs.

See run-sdp-amazon.jl for an example.

------------------------------------------------
LP Relaxation for O(log n) approximation
------------------------------------------------

run-LR-amazon.jl
run-LR-mathoverflow.jl
run-LR-trivago.jl

------------------------------
Hypergraph Cut Matching Games
------------------------------

run-HCM-(type).jl

where "type" is "amazon", "mathoverflow", or "trivago".

------------------------------
Plotting Figures
------------------------------

First need to run 

save-plotdata-mathoverflow.jl
save-plotdata-trivago.jl

and then plot using

plot-math-small.jl
plot-trivago-small.jl
