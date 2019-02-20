% Problem 4.(b), l1 ball constrained logistic regression
load("../datasets/rcv1.binary/rcv1_train.mat");

prox = @(x) ProjectOntoL1Ball(x, 10);
grad_g = @(x) (1/size(b, 1)*sum(-b.*A./(1+exp(-b.*(A*x))), 1))';

x0 = zeros(size(A, 2), 1);

%% For PG without line search

% L = norm(A'*A);
% L = 1.907300354573517e+05, eta <= 5.24e-6
eta = 1e-2;
epsilon = 1e-6;

xs = proximal_gradient(grad_g, prox, x0, eta, epsilon);

f = @(x) 1/size(b, 1)*sum(log(1+exp(-b.*(A*x))), 1);
fx = zeros(1, size(xs, 2));
for i = 1:size(xs, 2)
    fx(i) = f(xs(:, i));
end

semilogy(fx); hold on;

%% For PG with line search
eta = 1.0;
gamma_inc = 0.5;
gamma_dec = 0.5;
g = @(x) 1/2*norm(A*x-b)^2;

xs2 = proximal_gradient_with_line_search(g, grad_g, prox, x0, eta, gamma_inc, gamma_dec, epsilon);
fx2 = zeros(1, size(xs2, 2));
for i = 1:size(xs2, 2)
    fx2(i) = f(xs2(:, i));
end
semilogy(fx2);

legend('PG wo line search', "PG with line search");
