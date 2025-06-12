function [Wo,yo] = prox_tg(t,W,G,y,eps,rho)
  
    [k,l] = size(G);
    gradf_W = zeros(size(W));
    gradf_y = zeros(size(y));
    I = eye(size(W));
    for i = 1:k
        gradf_W = gradf_W + phiprim(sqrt(eps + norm(W*G(i,1:l-1)')^2) + G(i,1:l-1)*y - G(i,l))/(sqrt(eps + norm(W*G(i,1:l-1)')^2)).*((G(i,1:l-1)')*G(i,1:l-1));
        gradf_y = gradf_y + phiprim(sqrt(eps + norm(W*G(i,1:l-1)')^2) + G(i,1:l-1)*y - G(i,l))*G(i,1:l-1)';
    end
    gradf_W = W*gradf_W;
    W_b = W - t*gradf_W;
    [R,D] = eig((W_b - 2*eps*I + W_b')/2);

    W_ppp = zeros(size(D));
    for i = 1:length(D)
        if t + eps*rho*D(i,i) > 0
            W_ppp(i,i) = 0.5*(D(i,i) - eps + sqrt( (D(i,i) - eps)^2 + 4*(eps*D(i,i) + t/rho) ));
        else
            W_ppp(i,i) = 0;
        end
    end
    Wo = R'*W_ppp*R + eps*I;
    yo = y - t*gradf_y;
    
end