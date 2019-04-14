function xpath = level_set_with_extra_gradient(H, alpha, t0, x0, y0, epsilon, proj_x, proj_y, grad_x, grad_y, eta, K)
% Level set method with Extra Gradient being subprocess
% Input:
%   - H: function of (x, t), H(x, t) = max{f0(x)-t, f1(x), ..., fm(x)}.
%   - alpha: step size
%   - t0, x0, y0: initial guess
%   - epsilon: stop tolerance 
%   - proj_x, ..., K: parameters for extra gradient method

while(1)
    [x, y] = extra_gradient(proj_x, proj_y, grad_x, grad_y, x0, y0, eta, K);
    