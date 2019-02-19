% Problem 4.(a), box constrained linear regressiong
load("../datasets/E2006/E2006_train.mat");
grad_g = @(x) A'*(A*x-b);

l = -10*ones(size(A, 2), 1); u = 10*ones(size(A, 2), 1);
prox = @(x) proximal_box(x, l, u);

x0 = zeros(size(A, 2), 1);

% L = norm(A'*A);
% L = 1.907300354573517e+05, eta <= 5.24e-6
eta = 5e-6;
epsilon = 1e-6;

xs = proximal_gradient(grad_g, prox, x0, eta, epsilon);

f = @(x) 1/2*norm(A*x-b).^2;
fx = zeros(1, size(xs, 2));
for i = 1:size(xs, 2)
    fx(i) = f(xs(:, i));
end

semilogy(fx);


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