function [inform, x] = BFGS(fun, x, qnparams)
global numf numg;
numf = 0;
numg = 0;
lsparams = struct('c1', 1.0e-4, 'c2', 0.4, 'maxit', 20);

I = eye(size(x.p,1));
x.f = feval(fun, x.p, 1);
x.g = feval(fun, x.p, 2);
H = I;

for iter = 1 : qnparams.maxit
    p = -1 * H * x.g;
    [alpha, x_neo] = StepSize(fun, x, p, 1, lsparams);
    s = x_neo.p - x.p;
    y = x_neo.g - x.g;
    rho = 1.0 / (y' * s);
    if iter == 1
        H = (y' * s) / (y' * y) * I;
    end
    H = (I - rho * s * y') * H * (I - rho * y * s') + rho * s * s';
    x = x_neo;
    if norm(x.g) <= qnparams.toler * (1+abs(x.f))
        inform.status = 1;
        inform.iter = iter;
        return;
    end
end
inform.status = 0;
inform.iter = qnparams.maxit;
return;