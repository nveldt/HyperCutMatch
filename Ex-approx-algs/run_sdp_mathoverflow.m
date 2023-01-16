load ../data/mathoverflow-small/math_summary
diary Output/mathoverflow/Mathoverflow_SDPs_matlab_output.txt
for i = 4:numel(SetInds)
    num = SetInds(i);
    load(strcat('../data/mathoverflow-small/math_',num2str(num),'_lcc.mat'))
    fprintf('Mathoverflow dataset number %d \n\n',num)
    [X,Y, timeSDP, lowerbound,status] = hypergraphexpSDP(H,1800);
    fprintf('total time = %f \n',timeSDP)
    save(strcat('Output/mathoverflow/math_',num2str(num),'_SDP_solution.mat'),'X','Y','lowerbound','status','timeSDP')
end
diary off