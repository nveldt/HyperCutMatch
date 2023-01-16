load ../data/trivago-small/trivago_summary
diary Output/trivago/Trivago_SDPs_matlab_output.txt
for i = 9:numel(SetInds)
    num = SetInds(i);
    load(strcat('../data/trivago-small/trivago_',num2str(num),'_lcc.mat'))
    fprintf('Trivago dataset number %d \n\n',num)
    [X,Y, timeSDP, lowerbound,status] = hypergraphexpSDP(H,1800);
    fprintf('total time = %f \n',timeSDP)
    save(strcat('Output/trivago/trivago_',num2str(num),'_SDP_solution.mat'),'X','Y','lowerbound','timeSDP','status')
end
diary off