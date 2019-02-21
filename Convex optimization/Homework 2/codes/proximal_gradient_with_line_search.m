function [x] = proximal_gradient_with_line_search(g, grad_g, prox, x0, eta, gamma_inc, gamma_dec, epsilon)
% Proximal Gradient with line search; objective function f = g+h;
% Input:
%       g: a function of x
%       grad_g: the gradient of g
%       prox: the proximal mapping w.r.t eta*h
%       x0: initial guess
%       eta: initial step size
%       gamma_inc: increase rate of step size
%       gamma_dec: decrease rate of step size
%       epsilon: terminate condition

x = [x0];
while(1)
    eta = eta/gamma_inc;
    while(1)
        eta = eta*gamma_dec;
        
        x1 = prox(x0 - eta*grad_g(x0));
        x = [x x1];
        
        if(norm(x1 - x0) < epsilon)
            return;
        end
        
        size(x, 2)
        % more than 200 iterations
        if(size(x, 2) > 200)
            warning("May not converge");
            return;
        end
        
        if(g(x1) - g(x0) <= grad_g(x0)'*(x1-x0) + 1/(2*eta)*norm(x1-x0)^2)
            x0 = x1;
            break;
        end
        
        x0 = x1;
    end
end
