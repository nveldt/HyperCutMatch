function para = card_para(alpha, edgesize, type)
    if strcmp(type,'decreasing_rate')
        S = 1:edgesize/2;
        para = (1 - (2*(S-1)/(edgesize-1)).^alpha).^(1/alpha);
        if rem(edgesize, 2) == 0,
            para = [para -para(end:-1:1)];
        else
            para = [para 0 -para(end:-1:1)];
        end
    end
    if strcmp(type,'threshold_linear_rate'),
        threshold = ceil(alpha*edgesize);
        para = [ones(1, threshold) zeros(1, edgesize - 2*threshold) -ones(1, threshold)];
    end
    if strcmp(type,'threshold_quad_rate'),
        threshold = ceil(alpha*edgesize);
        S = 0:edgesize;
        para = S.*(edgesize-S)./(edgesize-1);
        para = min(para, threshold*(edgesize-threshold)/(edgesize-1));
        para = diff(para);
    end
    if strcmp(type,'threshold_linear_int'),
        para = [ones(1, alpha) zeros(1, edgesize - 2*alpha) -ones(1, alpha)];
    end
    if strcmp(type,'threshold_quad_int'),
        S = 0:edgesize;
        para = S.*(edgesize-S)./(edgesize-1);
        para = min(para, alpha*(edgesize-alpha)/(edgesize-1));
        para = diff(para);
    end

    if strcmp(type,'threshold_bound'),
        threshold = ceil(alpha(1)*edgesize);
        if threshold <= 1,
            para = zeros(1,edgesize);
            if edgesize > 1,
                para(1) = alpha(2);
                para(end) =  -alpha(2);
            end
            return;
        end
        S = 0:edgesize;
        para = S.*(edgesize-S)./(edgesize-1);
        bound = threshold*(edgesize-threshold)/(edgesize-1);
        para = min(para, bound);
        para = diff(para);
        if bound ~= 1, 
            para(2:end-1) = para(2:end-1)*(alpha(2)-1)/(bound-1);
        else
            para(2:end-1) = 0;
        end
    end
    if strcmp(type,'threshold_bound2'),
        threshold = ceil(alpha(1)*edgesize);
        S = 0:edgesize;
        para = S.*(edgesize-S)./(edgesize-1);
        bound = threshold*(edgesize-threshold)/(edgesize-1);
        para = min(para, bound);
        para = diff(para);
        para = para*alpha(2)/bound;
%         if bound ~= 1, 
%             para(2:end-1) = para(2:end-1)*(alpha(2)-1)/(bound-1);
%         else
%             para(2:end-1) = 0;
%         end
    end
end