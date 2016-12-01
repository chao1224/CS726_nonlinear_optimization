%% global variable
mu=0.01; L=1; kappa=L/mu;
n=100;
A=randn(n,n); [Q,R] = qr(A);
D=rand(n,1); D=10.^D; Dmin=min(D); Dmax=max(D);
D=(D-Dmin)/(Dmax-Dmin);
D=mu+D*(L-mu);
A=Q'*diag(D)*Q;

epoch_step  = 30000;
trial_num = 10;
av_sd_list=ones(trial_num,1)*epoch_step;
av_sde_list=ones(trial_num,1)*epoch_step;
av_nest_list=ones(trial_num,1)*epoch_step;
av_cg_list=ones(trial_num,1)*epoch_step;

for step = 1:trial_num
    x0=randn(n,1);%use a different x0 for each trial
   
    %% Steepest descent with sd \alpha 
    x_sd_list = zeros(n, epoch_step);
    x_sd_list(:,1)=x0;

    for k = 1:epoch_step-1
        alpha = 1.0 / L;
        x_k = x_sd_list(:,k);
        if 0.5 * x_k' * A * x_k <= 10^(-6)
            av_sd_list(step)=k;
            break
        end
        delta_k = A * x_k;
        x_sd_list(:,k+1) = x_k - alpha * delta_k;
    end


    %% Steepest descent with exact line search
    x_sde_list = zeros(n, epoch_step);
    x_sde_list(:,1)=x0;

    for k = 1:epoch_step-1
        x_k = x_sde_list(:,k);
        if 0.5 * x_k' * A * x_k <= 10^(-6)
            av_sde_list(step)=k;
            break
        end
        delta_k = A * x_k;
        alpha = (delta_k' * delta_k) / (delta_k' * A * delta_k) ;
        x_sde_list(:,k+1) = x_k - alpha * delta_k;
    end


    %% Nesterov's method
    x_nest_list = zeros(n, epoch_step);
    x_nest_list(:,1)=x0;

    m=mu;
    alpha=1.0/L;
    beta=(sqrt(L)-sqrt(m))/(sqrt(L)+sqrt(m));
    beta=(sqrt(kappa)-1)/(sqrt(kappa)+1);

    for k = 1:epoch_step-1
        x_k = x_nest_list(:,k);
        if 0.5 * x_k' * A * x_k <= 10^(-6)
            av_nest_list(step)=k;
            break
        end
        if k == 1
            y_k = x_k ;
        else
            y_k = x_k + beta * (x_k - x_nest_list(:,k-1));
        end

        delta_k = A * y_k;
        x_nest_list(:,k+1) = y_k - alpha * delta_k;
    end

    %% Conjugate gradient
    x_cg_list = zeros(n, epoch_step);
    x_cg_list(:,1)=x0;

    r = A * x0;
    p = -r;

    for k = 1:epoch_step-1
        x_k = x_cg_list(:,k);
        if 0.5 * x_k' * A * x_k < 10^(-6)
            av_cg_list(step)=k;
            break;
        end
        alpha = (r' * r) / (p' * A * p);
        x_cg_list(:,k+1) = x_k + alpha * p;
        r_prev = r;
        r = r + alpha * A * p;
        beta = r' * r / (r_prev' * r_prev);
        p = -r + beta * p;
    end
    
    
end

av_sd = mean(av_sd_list);
av_sde = mean(av_sde_list);
av_nest = mean(av_nest_list);
av_cg = mean(av_cg_list);
table = [[1:10]', av_sd_list, av_sde_list, av_nest_list, av_cg_list];

fprintf(1, ' steepest descend -fixed steps : %7.1f\n', av_sd);
fprintf(1, ' steepest descend -exact steps : %7.1f\n', av_sde);
fprintf(1, ' Nestrov : %7.1f\n', av_nest);
fprintf(1, ' conjugate gradient : %7.1f\n', av_cg);