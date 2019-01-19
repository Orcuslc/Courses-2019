function [x, v, t] = larmor_motion_solver(E, B, q, m, x0, v0, T0, T, dt)
% The 3D larmor motion equation is:
% (1) dx/dt = v
% (2) dv/dt = q/m*(E+ v \cross B)

% Input:
% E: electric field, a function of x and t
% B: magnectic field, a function of x and t
% q: charge
% m: mass
% x0, v0: initial conditions
% T0, T, dt: initial time, final time and timestep

% Output: 
% x, v: 3xL matrix, each column is a solution
% t: 1xL array, the timestamps

t = T0:dt:T;
x = zeros(3, length(t));
v = zeros(3, length(t));
x(:, 1) = x0; v(:, 1) = v0;

for i = 1:(length(t)-1)
    x(:, i+1) = x(:, i) + v(:, i)*dt;
    v(:, i+1) = v(:, i) + q/m*(E(x(:, i), t(i)) + cross(v(:, i), B(x(:, i), t(i))))*dt;
end
    