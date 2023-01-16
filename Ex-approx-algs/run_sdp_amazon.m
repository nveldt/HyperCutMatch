diary Output/amazon/Amazon_SDPs_matlab_output.txt
for num = 1:3
    load(strcat('../data/amazon-9/Am_',num2str(num),'_lcc.mat'))
    fprintf('Amazon dataset number %d \n\n',num)
    [X,Y, timeSDP, lowerbound,status] = hypergraphexpSDP(H,300);
    fprintf('total time = %f \n',timeSDP)
    save(strcat('Output/amazon/Am_',num2str(num),'_SDP_solution.mat'),'X','Y','lowerbound','timeSDP','status')
end
diary off