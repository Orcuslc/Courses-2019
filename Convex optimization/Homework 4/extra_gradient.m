function [xpath, ypath] = extra_gradient(proj_x, proj_y, grad_x, grad_y, x0, y0, eta, K)
% Extra gradient method for
% min_{x\in X} max_{y\in Y} f(x, y),
% where 
%   1. f is convex in x, concave in y, 
%   2. X, Y convex and compact
%   3. grad f is L-Lipschitz
% Input:
%   proj_x, proj_y: Projection operator onto X, Y
%   grad_x, grad_y: Partial derivative of f w.r.t. x, y
%   x0, y0: initial value
%   eta: step size, in (0, 1/L]
%   K: number of iterations

%% initialize
xpath = [x0];
ypath = [y0];
xsum = x0;
ysum = y0;

%% main iteration
for k = 1:K
    disp(k);
    u = proj_x(x0 - eta*grad_x(x0, y0));
    v = proj_y(y0 + eta*grad_y(x0, y0));
    x1 = proj_x(x0 - eta*grad_x(u, v));
    y1 = proj_y(y0 + eta*grad_y(u, v));
    
    xsum = xsum+x1;
    ysum = ysum+y1;
    xpath = [xpath xsum./(k+1)];
    ypath = [ypath ysum./(k+1)];
    
    x0 = x1;
    y0 = y1;
end