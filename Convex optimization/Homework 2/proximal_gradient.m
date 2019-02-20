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
    x1 = prox(x0 - eta*grad_g(x0));
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