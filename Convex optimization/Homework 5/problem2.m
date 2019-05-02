% Problem 2
% Augmented linearized Lagrangian method for Overlapped Group Regularized Logistic Regression
% Here f = g+h

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

% smoothing constant
mu = 0.001;

% gradient of g
y_star = @(x) ystar(x, C, lambda, mu);
grad_g = @(x) (1/size(b, 1)*sum(-M+M./(1+exp(-M*x)), 1))' + lambda*C'*y_star(x);

% proximal of h
prox = @(x, eta) x/max(10, norm(x, 2));

% initial guess
x0 = ones(54, 1)*0.01;

% step size 
% lipschitz constant
L = norm(C)/mu;
eta = 1/(2*L);
gamma0 = eta/2;
epsilon = 1e-6;
mu1 = 0;

x = accelerated_proximal_gradient(grad_g, prox, x0, eta, gamma0, mu1, epsilon);
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

% compute y_star for h_mu
function y = ystar(x, C, lambda, mu)

% projection operator onto unit L2 ball
S2 = @(x) x/max(1, norm(x, 2));

y = zeros(144, 1);
x = C*x;
for g = 1:16
    y(9*g-8:9*g) = S2(lambda/mu*x(9*g-8:9*g));
end
end

function y = target(x, M, lambda)
y = 1/size(M, 1)*sum(log(1+exp(-M*x)), 1);
for g = 1:16
    y = y + lambda*norm(x(3*g-2:3*g+6), 2);
end
end
