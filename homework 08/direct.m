function [inform, x] = direct(fun, x, directparams)
global numf numg;
numf = 0;
numg = 0;
n = size(x.p, 1);
directions = [eye(n), -1*eye(n)];
gamma = 1;
x.f = feval(fun, x.p, 1);

for step = 1:directparams.maxit
    order = randperm(2*n);
    flag = 0;
    for k = 1 : 2*n
        d = directions(:, order(k));
        x_neo = struct('p', x.p + gamma*d);
        x_neo.f = feval(fun, x_neo.p, 1);
        if x_neo.f < x.f - gamma^2
            x = x_neo;
            gamma = gamma * directparams.phi;
            flag = 1;
            break;
        end
    end
    if flag == 0
        gamma = gamma * directparams.theta;
    end
    if gamma <= directparams.toler
        inform.status = 1;
        inform.iter = step;
        return;
    end
end

inform.status = 0;
inform.iter = directparams.maxit;
return