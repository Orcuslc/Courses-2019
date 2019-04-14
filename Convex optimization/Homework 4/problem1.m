% Problem 1, 
% Extra gradient for overlapping group regularized logistic regression

load("../datasets/covtype/covtype.mat");

% convert b to {-1, 1}
b = 2*b-3;

% Group matrix
C = zeros(144, 54);
for g = 1:16
    C(9*g-8:9*g, 3*g-2:3*g+6) = eye(9);
end


% compute b.*A
M = diag(sparse(b))*A;

% penalty 
lambda = 0.0001;

% partial derivative
grad_x = @(x, y) (1/size(b, 1)*sum(-M+M./(1+exp(-M*x)), 1))' + lambda*C'*y;
grad_y = @(x, y) lambda*C*x;

% projection operator for X and Y
proj_x = @(x) x/max(10, norm(x, 2));
proj_y = @Proj_Y;

% step size
eta = 0.00025;

% initial guess
x0 = ones(54, 1)*0.1;
y0 = ones(144, 1)*0.1;

% iteration
K = 200;
[x, y] = extra_gradient(proj_x, proj_y, grad_x, grad_y, x0, y0, eta, K);

% target function
f = @(x) target(x, M, lambda);
fx = zeros(1, size(x, 2));
for i = 1:size(x, 2)
    fx(i) = f(x(:, i));
end
semilogy(fx);
xlabel("K");
ylabel("function value");
title("objective value w.r.t. average iterate");

% projection operator for Y
function y = Proj_Y(y)
for g = 1:16
    y(9*g-8:9*g) = y(9*g-8:9*g)/max(1, norm(y(9*g-8:9*g), 2));
end
end

function y = target(x, M, lambda)
y = 1/size(M, 1)*sum(log(1+exp(-M*x)), 1);
for g = 1:16
    y = y + lambda*norm(x(3*g-2:3*g+6), 2);
end
end