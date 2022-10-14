function [labels, incidence_partition, NCut] = optthreshold(incidence_list, parameter_list, submodular_type, mu, N, R, x)
    [xsort, index] = sort(x, 'descend');
    xtemp = zeros(1, length(index));
    sum_mu = sum(mu);
    NCut = inf;
    optpos = 0;
    for j = 1:length(index) - 1,
        xtemp(index(j)) = 1; 
        fx = 0;
        for i = 1:length(incidence_list),
            if strcmp(submodular_type{i},'edge'),
                fx = fx + parameter_list{i}(1) * abs(xtemp(incidence_list{i}(1))-xtemp(incidence_list{i}(2)));
            end
            if strcmp(submodular_type{i},'para_edge'),
                templist = reshape(incidence_list{i}, 2, length(parameter_list{i}));
                fx = fx + sum(parameter_list{i} .* abs(xtemp(templist(1,:))-xtemp(templist(2,:))));
            end
            if strcmp(submodular_type{i},'hyperedge'),
                tempx = xtemp(incidence_list{i});
                fx = fx + parameter_list{i}(1) * (max(tempx)-min(tempx));
            end
            if strcmp(submodular_type{i},'concave_card'),
                cardS = sum(xtemp(incidence_list{i}));
                fx = fx + sum(parameter_list{i}(1:cardS));
            end
        end
        vol = sum(xtemp.*mu);
        fx = fx/min(vol, sum_mu - vol);
        if fx < NCut,
            NCut = fx;
            optpos = j;
        end
    end
    labels = zeros(1, length(index));
    labels(index(1:optpos)) = 1;
    incidence_partition = zeros(1, length(incidence_list));
    for i = 1:length(incidence_list),
        templabels = labels(incidence_list{i});
        if max(templabels) ~= min(templabels),
            incidence_partition(i) = 2;
        else
            incidence_partition(i) = templabels(1);
        end
    end    
end