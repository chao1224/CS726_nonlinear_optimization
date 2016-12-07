x = [0 0; 1 0; 2 0; 1 1; 0 2; 0 1];
y = [1; 2.0084; 7.0091; 1.0168; -0.9909; -0.9916];


%% initial
n = size(x, 2);
N = (n+1)*(n+2) / 2;
big = zeros(N, N);

%% transform
for row = 1 : N
    big(row, 1) = 1;
    big(row, 2:1+n) = x(row,:);
    col = 1+n;
    for i = 1:n
        for j = i:n
            col = col + 1;
            if i == j
                big(row, col) = 0.5*x(row,i)^2;
            else
                big(row, col) = x(row,i) * x(row, j);
            end
        end
    end
end

%% calculate

solution = pinv(big) * y;

%% transform back

c = solution(1);
g = zeros(n,1);
g(1:n) = solution(2:1+n);
G = zeros(n,n);
col = n+1;
for i = 1:n
    for j = i:n
        col = col + 1;
        G(i, j) = solution(col);
        G(j, i) = solution(col);
    end
end