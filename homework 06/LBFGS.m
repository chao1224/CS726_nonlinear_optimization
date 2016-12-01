 function [inform, x] = LBFGS(fun, x, lbfgsparams)
global numf numg;
numf = 0;
numg = 0;
lsparams = struct('c1', 1.0e-4, 'c2', 0.4, 'maxit', 20);

n = size(x.p,1);
m = lbfgsparams.m;
I = eye(n);
x.f = feval(fun, x.p, 1);
x.g = feval(fun, x.p, 2);
y = 1;
s = 1;

s_list = zeros(n, m);
y_list = zeros(n, m);
rho_list = zeros(1, m);
alpha_list = zeros(1, m);

for iter = 1 : lbfgsparams.maxit
    %% calculate gradient direction
    H = (s' * y) / (y' * y) * I;
    q = x.g;
    for i = iter-1: -1 :max(iter-m,1)
        s_t = s_list(:,mod(i,m)+1);
        y_t = y_list(:,mod(i,m)+1);
        rho_t = rho_list(mod(i,m)+1);
        alpha = rho_t * s_t' * q;
        alpha_list(mod(i,m)+1) = alpha;
        q = q - alpha * y_t;
    end
    r = H * q;
    for i = max(iter-m,1) : iter-1
        s_t = s_list(:,mod(i,m)+1);
        y_t = y_list(:,mod(i,m)+1);
        rho_t = rho_list(mod(i,m)+1);
        beta = rho_t * y_t' * r;
        alpha = alpha_list(mod(i,m)+1);
        r = r + s_t * (alpha - beta);
    end
    p = -1 * r;
    %% calculate next x
    [alpha, x_neo] = StepSize(fun, x, p, 1, lsparams);
    %% update iter-m
    s = x_neo.p - x.p;
    y = x_neo.g - x.g;
    rho = 1.0 / (y' * s);
    s_list(:, mod(iter,m)+1) = s;
    y_list(:, mod(iter,m)+1) = y;
    rho_list(mod(iter,m)+1) = rho;
    x = x_neo;
    %% decide if terminate
    if norm(x.g) <= lbfgsparams.toler * (1+abs(x.f))
        inform.status = 1;
        inform.iter = iter;
        return;
    end
end
inform.status = 0;
inform.iter = lbfgsparams.maxit;
return;