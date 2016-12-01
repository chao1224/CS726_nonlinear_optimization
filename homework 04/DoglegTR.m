 function [inform, x] = DoglegTR(fun, x, trparams)
numf = 0;
numg = 0;
numH = 0;
delta_cur = trparams.Delta0;

x.f = feval(fun, x.p, 1);
x.g = feval(fun, x.p, 2);
if norm(x.g) <= trparams.toler
    inform.status = 1;
    inform.iter = 0;
    return
end
x.h = feval(fun, x.p, 4);
[v,d] = eig(x.h);
for i = 1 : size(d,1)
    d(i,i) = max(d(i,i),trparams.delta);
end
x.h = v' * d * v;

for ite = 1 : trparams.maxit
    inform.iter = ite;
    %% update with dogleg
    p_B = - x.h^(-1) * x.g;
    p_u = - (x.g' * x.g * x.g )  / (x.g' * x.h * x.g);
    if norm(p_B) <= delta_cur
        p_hat = p_B;
    elseif norm(p_u) >= delta_cur
        p_hat = p_u * delta_cur / norm(p_u);
    else
        a = 1.0 * norm(p_B-p_u).^2;
        b = 1.0* 2 *  p_u' * (p_B -p_u);
        c = 1.0 * norm(p_u).^2 - delta_cur.^2;
        tao1 =  (-b + sqrt(b*b - 4*a*c)) / (2*a);
        tao2 =  (-b - sqrt(b*b - 4*a*c)) / (2*a);
        if 0 <= tao1 && tao1 <= 1
            tao = tao1;
        else
            tao = tao2;
        end
        p_hat = p_u + (tao) * (p_B - p_u);
    end
    p_neo = p_hat;
    
    next_f = feval(fun, x.p+p_neo,1);
    rho_cur = 1.0 * (x.f - next_f) / -(x.g' * p_neo + 0.5 * p_neo' * x.h * p_neo );
    if rho_cur < 0.25
        delta_next = 0.25 * delta_cur;
    elseif rho_cur > 0.75 && norm(p_neo) == delta_cur
        delta_next = min(2*delta_cur,  trparams.hatDelta);
    else
        delta_next = delta_cur;
    end
    delta_cur = delta_next;
    
    if rho_cur > trparams.eta
        x.p = x.p + p_neo;
        x.f = next_f;
        x.g = feval(fun, x.p, 2);
       %% check if stopping criteria meets
        if norm(x.g) <= trparams.toler
            inform.status = 1;
            return
        end
        x.h = feval(fun, x.p, 4);
        [v,d] = eig(x.h);
        for i = 1 : size(d,1)
            d(i,i) = max(d(i,i),trparams.delta);
        end
        x.h = v' * d * v;
    end
    fprintf(1, ' iter %3d: f=%12.5e, ||Df||=%12.5e, Delta=%7.2e\n', inform.iter, x.f, norm(x.g), delta_cur);

    
end

inform.status = 0;
inform.iter = trparams.maxit;

return;