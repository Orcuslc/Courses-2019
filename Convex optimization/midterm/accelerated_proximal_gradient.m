function [x] = accelerated_proximal_gradient(grad_g, prox, x0, eta, gamma0, mu, epsilon)
% Accelerated Proximal Gradient without line search; objective function f = g+h;
% Input:
%       grad_g: the gradient of g
%       prox: the proximal mapping w.r.t eta*h;
%       x0: initial guess
%       eta: step size
%       gamma0: initial control variable
%       mu: the strong convex parameter
%       epsilon: terminate condition

x = [x0];
z = x0;

while(1)
    alpha = (eta/2)*(-(gamma0-mu)+sqrt((gamma0-mu)^2+4*gamma0/eta));
    gamma1 = alpha^2/eta;
    y = 1/(alpha*gamma0 + gamma1)*(alpha*gamma0*z + gamma1*x0);
    x1 = prox(y - eta*grad_g(y), eta);
    disp(norm(grad_g(y)));
    z = x0 + 1/alpha*(x1 - x0);
    
    x = [x x1];
    if(norm(x1 - x0) < epsilon)
        return;
    end
    x0 = x1;
    
    disp(size(x, 2))
    % more than 200 iterations
    if(size(x, 2) > 200)
        warning("May not converge");
        return;
    end
end