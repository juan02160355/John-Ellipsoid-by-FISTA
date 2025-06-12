%% Given
constraint = [0 -0.25 1;-0.5*sqrt(3) 0.5 1;0.5*sqrt(3) 0.5 1;-0.5*sqrt(3) -0.5 0.8;0 1 1];%Edit constraint here
eps = 10^-16;
rho = 40;
beta = 0.401;
t = 19;
[K,l] = size(constraint);

%% Initialize
k = 1;
u_1 = 0;
u_2 = 0;
W = [0 0;0 0];
y = [0;0];
W_1 = [1 0;0 1];
W_2 = [1 0;0 1];
Wk = [0 0;0 0];
y_1 = [1;0];
y_2 = [0;1];
yk = [0;0];


%% FISTA
while true
    Wk_origin = Wk;
    yk_origin = yk;

    u_1 = 0.5*(1+sqrt(1+4*u_2^2));
    W = W_1 + ((u_2-1)/u_1)*(W_1-W_2);
    y = y_1 + ((u_2-1)/u_1)*(y_1-y_2);
    [Wk,yk] = prox_tg(t,W,constraint,y,eps,rho);

    gradf_W = zeros(size(W));
    gradf_y = zeros(size(y));
    
    for i = 1:K
        gradf_W = gradf_W + phiprim(sqrt(eps + norm(W*constraint(i,1:l-1)')^2) + constraint(i,1:l-1)*y - constraint(i,l))/(sqrt(eps + norm(W*constraint(i,1:l-1)')^2)).*((constraint(i,1:l-1)')*constraint(i,1:l-1));
        gradf_y = gradf_y + phiprim(sqrt(eps + norm(W*constraint(i,1:l-1)')^2) + constraint(i,1:l-1)*y - constraint(i,l))*constraint(i,1:l-1)';
    end
    gradf_W = W*gradf_W;
    
    while diff_f(Wk,yk,constraint,eps) > diff_f(W,y,constraint,eps) + MIP([gradf_W,gradf_y],[Wk-W,yk-y]) + 1/(2*t)*norm([Wk-W,yk-y],'fro')^2
        t = beta*t;
        [Wk,yk] = prox_tg(t,W,constraint,y,eps,rho);
    end
    k = k + 1;
    
    W_2 = W_1;
    y_2 = y_1;
    W_1 = Wk;
    y_1 = yk;

    %stopping criteria
    if norm([Wk-Wk_origin,yk-yk_origin],'fro') < 0.0001
        break;
    end
end


%% Plot
n = 1:360;
n = (n/180)*pi;
x = cos(n);
u = sin(n);
a = [x;u];
e = zeros(2,360);
for i = 1:360
    e(:,i) = Wk*a(:,i) + y;
end

n = 1:101;
x = -5 + (n - 1)/10;

plot(e(1,:),e(2,:))
hold("on")
cy = zeros(K,101);
for i = 1:K
    cy(i,:) = (constraint(i,l) - x*constraint(i,1))/constraint(i,2);
    plot(x,cy(i,:));
end




