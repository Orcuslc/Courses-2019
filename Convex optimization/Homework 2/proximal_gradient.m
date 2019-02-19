function [x] = proximal_gradient(grad_g, prox, x0, eta, epsilon)
% Proximal Gradient without line search; objective function f = g+h;
% Input:
%       grad_g: the gradient of g
%       prox: the proximal mapping w.r.t eta*h
%       x0: initial guess
%       eta: step size
%       epsilon: terminate condition

x = [x0];
while(1)
    x1 = prox(x(end) - eta*grad_g(x(end)));
    x = [x, x1];
    if(norm(x1 - x(end-1)) < epsilon)
        return;
    end
end