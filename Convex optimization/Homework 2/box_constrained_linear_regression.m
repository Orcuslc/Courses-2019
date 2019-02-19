% Problem 4.(a), box constrained linear regression
A = ;
b = ;
grad_g = @(x) A*x-b;

l = -10*ones(size(b)); u = 10*ones(size(b));
prox = @(x) proximal_box(x, l, u);

x0 = zeros(size(b));
eta = 1.0;
epsilon = 1e-6;

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