function para = alpha_linear(alpha, edgesize)
    S = edgesize;
    f = zeros(S+1,1);
    for k = 0:S
        a = k/(ceil(alpha*S));
        b = (S-k)/(ceil(alpha*S));
        f(k+1) = 1/2 + 1/2*min([1,a,b]);
    end
    f(1) = 0;
    f(end) = 0;

    v = zeros(S,1);
    for j = 1:S
        v(j) = f(j+1) - f(j);
    end
    para = v';
end