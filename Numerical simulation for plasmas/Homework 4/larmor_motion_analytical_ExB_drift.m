function [x, v, t] = larmor_motion_analytical_ExB_drift(N)
% Analytical solution of larmor motion, with ExB drift of uniform magnetic
% and perpendicular electric field
% equation:
%   dx/dt = v
%   m(dv/dt) = q(E + v \cross B)

% N: number of timesteps

% Electromagnetic Field Input
B0 = [0. 0. 1.]';    % Initial Magnetic Field
E0 = [0. 0.1 0.]';   % Initial Electric Field
Bmag = norm(B0);     % Magnitude of B0

% Parameters of the ion
q = 1.;  % Charge
m = 1.;  % Mass

% Time span 
T0 = 0.;
T = 20*pi;
dt = (T-T0)/N;
t = (T0:dt:T)';

% Initial Conditions
x0 = [0. 1. 0. ]';  % Initial position
v0 = [1. 0. 0. ]';  % Initial velocity

% Analytical Solution
x = zeros(3, size(t, 1));
v = zeros(3, size(t, 1));

% seperate the velocity into two parts: v = ve + u, where ve is the
% constant solution of the equation, i.e., making RHS = 0.
drift = cross(E0, B0)/Bmag^2;
u = v0 - drift;
% u = v0;
u_perpendicular = sqrt(u(1)^2+u(2)^2);
u_parallel = u(2);
omega = q*B0(3)/m;

v(1, :) = u_perpendicular*cos(omega*t);
v(2, :) = -u_perpendicular*sin(omega*t);
v(3, :) = q/m*E0(3)*t+u(3);
v = v + drift;

x(1, :) = u_perpendicular/omega*sin(omega*t) + drift(1)*t + x0(1);
x(2, :) = u_perpendicular/omega*cos(omega*t) + drift(2)*t + x0(2)-u_perpendicular/omega;
x(3, :) = 0.5*q/m*E0(3)*(t.^2) + u(3)*t + x0(3);