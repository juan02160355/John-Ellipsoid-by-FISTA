function out = diff_f(W,y,G,eps)
    out = 0;
    [K,l] = size(G);
    for i = 1:K
        out = out + Hb( sqrt(eps + norm(W*G(i,1:l-1)')^2) + G(i,1:l-1)*y - G(i,l));
    end
end