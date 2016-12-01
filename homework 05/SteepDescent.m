function [inform, x] = SteepDescent(fun, x, sdparams)

stepSizeParam = struct('c1',0.01, 'c2', 0.3, 'maxit',100);
stepSizeParam = struct('c1',0.001, 'c2', 0.5, 'maxit',100);
x.f = feval(fun, x.p, 1);
x.g = feval(fun, x.p, 2);

for i = 1 : sdparams.maxit
    [alpha, x] = StepSize(fun, x, -x.g, 1, stepSizeParam);
    if norm(x.g, Inf) <= sdparams.toler * (1+abs(x.f))
        inform.status = 1;
        inform.iter = iter;
        return;
    end
end

inform.status = 0;
inform.iter = sdparams.maxit;
return;