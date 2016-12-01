n=500;
x_axis = [1:n];

y_sd = zeros(1,n);
y_sde = zeros(1,n);
y_nest = zeros(1, n);
y_cg = zeros(1,n);

for i = 1 : n
    y_sd(i) = log10(0.5 * x_sd_list(:,i)' * A * x_sd_list(:,i));
    y_sde(i) = log10(0.5 * x_sde_list(:,i)' * A * x_sde_list(:,i));
    y_nest(i) = log10(0.5 * x_nest_list(:,i)' * A * x_nest_list(:,i));
    y_cg(i) = log10(0.5 * x_cg_list(:,i)' * A * x_cg_list(:,i));
end

%ylim([10^(-6),2]);
plot(x_axis, y_sd, 'r');
hold on;
plot(x_axis, y_sde, 'g');
hold on;
plot(x_axis, y_nest, 'b');
hold on;
plot(x_axis, y_cg, 'm');
hold on;
title('Convergence Rate');
xlabel('iteration');
ylabel('log_{10}(f(x)-f(x^*))');
legend('SD:const', 'SD:exact', 'Nesterov', 'CG');
