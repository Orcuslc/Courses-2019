function [x, y] = primal_dual_subgradient(gx, gy, proj_x, proj_y, eta, x0, y0, K)
% Primal-Dual subgradient method for min_x max_y f(x, y)
% gx, gy: partial_x f, partial_y f
% proj_x, proj_y: projection operator w.r.t. X, Y
% eta: step size, a function of k
% x0, y0: initial guess
% K: number of iteration 

x = [x0];
y = [y0];

xsum = x0;
ysum = y0;

for k = 1:K
    disp(["iteration No.", k]);
    eta_k = eta(k);
    x1 = proj_x(x0 - eta_k*gx(x0, y0));
    y1 = proj_y(y0 + eta_k*gy(x0, y0));
    
    xsum = xsum + x1;
    ysum = ysum + y1;
    
    x = [x xsum/(k+1)];
    y = [y ysum/(k+1)];
    
    x0 = x1;
    y0 = y1;
end