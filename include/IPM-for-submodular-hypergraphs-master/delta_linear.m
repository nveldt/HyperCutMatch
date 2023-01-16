function para = delta_linear(delta, edgesize)
    S = edgesize;
    f = zeros(S+1,1);
    for k = 0:S
        f(k+1) = min(min(k,S-k),delta);
    end
    f(1) = 0;
    f(end) = 0;

    v = zeros(S,1);
    for j = 1:S
        v(j) = f(j+1) - f(j);
    end
    para = v';
end