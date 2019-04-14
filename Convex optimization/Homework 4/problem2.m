% Problem 2
% Two-person-zero-sum
% by Extra Gradient and Primal-Dual subgradient

load("two-person-zero-sum.mat");

proj_x = @(x) projection_on_simplex(x);
proj_y = @(y) projection_on_simplex(y);

grad_x = @(x, y) A'*y;
grad_y = @(x, y) A*x;

x0 = zeros(size(A, 2), 1); x0(1) = 1.0;
y0 = zeros(size(A, 1), 1); y0(1) = 1.0;

eta = 0.001;
K = 200;

%% Extra Gradient 
[x, y] = extra_gradient(proj_x, proj_y, grad_x, grad_y, x0, y0, eta, K);

% target function
f = @(x, y) target(A, x, y);
fx = zeros(1, size(x, 2));
for i = 1:size(x, 2)
    fx(i) = f(x(:, i), y(:, i));
end
figure;
plot(fx);
xlabel("K");
ylabel("function value");
title("Extra Gradient, objective value w.r.t. average iterate");


%% Primal-Dual subgradient
c = 0.005;
eta = @(k) c/sqrt(k);
[x, y] = primal_dual_subgradient(grad_x, grad_y, proj_x, proj_y, eta, x0, y0, K);
fx = zeros(1, size(x, 2));
for i = 1:size(x, 2)
    fx(i) = f(x(:, i), y(:, i));
end
figure;
plot(fx);
xlabel("K");
ylabel("function value");
title("Primal-Dual, objective value w.r.t. average iterate");

function v = target(A, x, y)
v = y'*A*x;
end
