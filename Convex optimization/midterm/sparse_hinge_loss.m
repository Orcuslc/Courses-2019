load('../datasets/news20/news20.mat');


c = -b.*A;
lambda = 0.01;
C = 0.1;
N = 19996;
px = @(x, eta) proximal_l1_norm(x, eta, lambda);
g = @(x, n) G(x, n, c);
x0 = zeros(size(A, 2), 1);
xs = stochastic_subgradient(g, px, C, N, x0);

f = @(x) mean(max(1+c*x, 0), 1) + lambda*norm(x, 1);
fs = zeros(1, size(xs, 2));
for i = 1:size(xs, 2)
    fs(i) = f(xs(:, i));
end
semilogy(fs);

function dx = G(x, n, c)
    z = -c(n, :);
    y = 1-z*x;
    if(y == 0)
        dx = rand();
    elseif(y > 0)
        dx = 1;
    else
        dx = 0;
    end
    dx = z'*dx;
end

function px = proximal_l1_norm(x, eta, lambda)

px = sign(x).*max(0, abs(x)-eta*lambda);

end