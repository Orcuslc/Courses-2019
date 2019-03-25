% Problem 4.(a), box constrained linear regression
load("../datasets/E2006/E2006_train.mat");
grad_g = @(x) A'*(A*x-b);

l = -10*ones(size(A, 2), 1); u = 10*ones(size(A, 2), 1);
prox = @(x) proximal_box(x, l, u);

x0 = zeros(size(A, 2), 1);

%% For PG without line search

% L = norm(A'*A);
% L = 1.907300354573517e+05, eta <= 5.24e-6
eta = 5e-6;
epsilon = 1e-6;

xs = proximal_gradient(grad_g, prox, x0, eta, 1/(2*eta), 0, epsilon);

f = @(x) 1/2*norm(A*x-b).^2;
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

legend('PG without line search', "PG with line search");

function px = proximal_box(x, l, u)
% proximal of a h(x) = 1_c(x), where c is a box constrained by l and u.
px = x;
for i = 1:length(x)
    if(x(i) > u(i))
        px(i) = u(i);
    elseif(x(i) < l(i))
        px(i) = l(i);
    end
end
end