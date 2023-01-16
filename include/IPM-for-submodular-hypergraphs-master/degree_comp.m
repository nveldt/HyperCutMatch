function degvec = degree_comp(incidence_list, N, type)
    degvec = zeros(1, N);
    for i = 1:length(incidence_list),
        degvec(incidence_list{i}) = degvec(incidence_list{i}) + 1;
    end
    if strcmp(type, 'partial'),
        new_deg_vec = cell(length(incidence_list),1);
        for i = 1:length(incidence_list),
            new_deg_vec{i} = degvec(incidence_list{i});
        end
        degvec = new_deg_vec;
    end
end