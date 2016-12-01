function [inform, x] = CG_FR(fun, x, nonCGparams)

stepSizeParam = struct('c1',0.01, 'c2', 0.3, 'maxit',100);
stepSizeParam = struct('c1',0.001, 'c2', 0.5, 'maxit',100);
x.f = feval(fun, x.p, 1);
x.g = feval(fun, x.p, 2);
p = -x.g;

for iter = 1 : nonCGparams.maxit
    [alpha, x_neo] = StepSize(fun, x, p, 1, stepSizeParam);
    beta = x_neo.g' * x_neo.g / (x.g' * x.g);
    p = -x_neo.g + beta * p;
    x = x_neo;
    if norm(x.g, Inf) <= nonCGparams.toler * (1+abs(x.f))
        inform.status = 1;
        inform.iter = iter;
        return;
    end
end

inform.status = 0;
inform.iter = nonCGparams.maxit;
return;