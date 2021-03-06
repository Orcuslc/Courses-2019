function [x, v, t] = larmor_motion_dimensionless_solver(E, B, x0, v0, T0, T, dt, order)
% The 3D dimensionless larmor motion equation is:
% (1) dx'/dt' = v'
% (2) dv'/dt' = E'+V'\cross B',
% where the normalizations are:
% x' = x/r_L,
% t' = Omega*t
% v' = v/vp (vp: perpendicular part)
% B' = B/B0
% E' = E/(vp*B0)
% and the angular ion cyclotron frequency and Larmor radius given by:
% Omega = q*B0/m,
% r_L = vp/Omega

% Input: (ALL INPUT ARE NORMALIZED)
% E: electric field, a function of x and t
% B: magnectic field, a function of x and t
% x0, v0: initial conditions
% T0, T, dt: initial time, final time and timestep
% order: 1 (Euler) or 2 (leapfrog)

% Output: 
% x, v: 3xL matrix, each column is a solution
% t: 1xL array, the timestamps

t = T0:dt:T;
x = zeros(3, length(t));
v = zeros(3, length(t));
x(:, 1) = x0; v(:, 1) = v0;

if nargin == 7
    order = 1; % default: Euler
end

switch order
    case 1 
        % Euler
        for i = 1:(length(t)-1)
            x(:, i+1) = x(:, i) + v(:, i)*dt;
            v(:, i+1) = v(:, i) + (E(x(:, i), t(i)) + cross(v(:, i), B(x(:, i), t(i))))*dt;
        end
    case 2
        % Leapfrog
        % the first step is to use Euler method and get x1 from x0, then
        % use leapfrog in iteration.
        x(:, 2) = x(:, 1) + v(:, 1)*dt;
        v(:, 2) = v(:, 1) + (E(x(:, 1), t(1)) + cross(v(:, 1), B(x(:, 1), t(1))))*dt;
        for i = 2:(length(t)-1)
           x(:, i+1) = x(:, i-1) + 2*v(:, i)*dt;
           v(:, i+1) = v(:, i-1) + 2*(E(x(:, i), t(i)) + cross(v(:, i), B(x(:, i), t(i))))*dt;
        end
end
        