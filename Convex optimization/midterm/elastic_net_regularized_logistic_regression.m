% Problem 4.
load('../datasets/news20/news20.mat');

% In order to guarantee the strong convexity, we need to regard this
% problem as:
% g(x) = logistic + \lambda_2/2*|x|_2^2
% h(x) = \lambda_1*|x|_1

% L = 1/2*lambda_max(A'A) = 585;
eta = 0.0015; 

lambda1 = 0.1;
lambda2 = 0.001;
c = -b.*A;
grad_g = @(x) (1/size(b, 1)*sum(c./(1+exp(-b.*(A*x))), 1))' + lambda2*x;
prox = @(x, eta) proximal_l1_norm(x, eta, lambda1);

% mu = lambda_2
mu = 0.001;
gamma0 = 2.0;
epsilon = 1e-6;

x0 = zeros(size(A, 2), 1);

xs = accelerated_proximal_gradient(grad_g, prox, x0, eta, gamma0, mu, epsilon);
f = @(x) 1/size(b, 1)*sum(log(1+exp(-b.*(A*x))), 1) + lambda1*norm(x, 1) + lambda2/2*norm(x, 2)^2;
fx = zeros(1, size(xs, 2));
for i = 1:size(xs, 2)
    fx(i) = f(xs(:, i));
end

semilogy(fx); hold on;



% function px = proximal_elastic_net(x, eta, lambda1, lambda2)
% % Get the form in 
% % https://stats.stackexchange.com/questions/236753/shrinkage-operator-for-elastic-net-regularization
% % and Page 7 of
% % https://www.cs.cmu.edu/~suvrit/teach/yaoliang_proximity.pdf
% 
% px = 1/(1+eta*lambda2)*proximal_l1_norm(x, eta, lambda1);
% 
% end

function px = proximal_l1_norm(x, eta, lambda)

px = sign(x).*max(0, abs(x)-eta*lambda);

end