%% global variable
mu=0.01; L=1; kappa=L/mu;
n=100;
A=randn(n,n); [Q,R] = qr(A);
D=rand(n,1); D=10.^D; Dmin=min(D); Dmax=max(D);
D=(D-Dmin)/(Dmax-Dmin);
D=mu+D*(L-mu);
A=Q'*diag(D)*Q;

epoch_step  = 30000;
x0=randn(n,1);%use a different x0 for each trial

%% Steepest descent with sd \alpha 
x_sd_list = ones(n, epoch_step);
x_sd_list(:,1)=x0;
av_sd=epoch_step;

for k = 1:epoch_step-1
    alpha = 1.0 / L;
    x_k = x_sd_list(:,k);
    if 0.5 * x_k' * A * x_k <= 10^(-6)
        av_sd=k;
        break
    end
    delta_k = A * x_k;
    x_sd_list(:,k+1) = x_k - alpha * delta_k;
end


%% Steepest descent with exact line search
x_sde_list = ones(n, epoch_step);
x_sde_list(:,1)=x0;
av_sde = epoch_step;

for k = 1:epoch_step-1
    x_k = x_sde_list(:,k);
    if 0.5 * x_k' * A * x_k <= 10^(-6)
        av_sde=k;
        break
    end
    delta_k = A * x_k;
    alpha = (delta_k' * delta_k) / (delta_k' * A * delta_k) ;
    x_sde_list(:,k+1) = x_k - alpha * delta_k;
end


%% Nesterov's method
x_nest_list = ones(n, epoch_step);
x_nest_list(:,1)=x0;
av_nest=epoch_step;

m=mu;
alpha=1.0/L;
beta=(sqrt(L)-sqrt(m))/(sqrt(L)+sqrt(m));

for k = 1:epoch_step-1
    x_k = x_nest_list(:,k);
    if 0.5 * x_k' * A * x_k <= 10^(-6)
        av_nest=k;
        break
    end
    if k == 1
        y_k = x_k + beta * x0;
    else
        y_k = x_k + beta * (x_k - x_nest_list(:,k-1));
    end
    
    delta_k = A * y_k;
    x_nest_list(:,k+1) = y_k - alpha * delta_k;
end

%% Conjugate gradient
x_cg_list = ones(n, epoch_step);
x_cg_list(:,1)=x0;
av_cg=epoch_step;
r = A * x0;
p = -r;

for k = 1:epoch_step-1
    x_k = x_cg_list(:,k);
    if 0.5 * x_k' * A * x_k < 10^(-6)
        av_cg=k;
        break;
    end
    alpha = - (r' * p) / (p' * A * p);
    x_cg_list(:,k+1) = x_k + alpha * p;
    r = A * x_k;
    beta = r' * A * p / (p' * alpha * p);
    p = -r + beta * p;
end
