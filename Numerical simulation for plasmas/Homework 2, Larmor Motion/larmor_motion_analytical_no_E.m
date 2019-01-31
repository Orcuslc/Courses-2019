function [x, v, t] = larmor_motion_analytical_no_E(N)
% Analytical solution of larmor motion, in Constant, Uniform B with E = 0.
% N: number of timesteps
% equation:
%   dx/dt = v
%   m(dv/dt) = qv \cross B

% Electromagnetic Field Input
B0 = [0. 0. 1.]';   % Initial Magnetic Field
E0 = [0. 0. 0.]';   % Initial Electric Field
Bmag = norm(B0);    % Magnitude of B0

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
v_perpendicular = sqrt(v0(1)^2+v0(2)^2);  % perpendicular component of velocity (perpendicular to B)
v_parallel = v0(3);                       % parallel component of velocity
rl = m*v_perpendicular/(q*B0(3));         % Larmor radius
omega = q*B0(3)/m;                        % Cyclotron frequency

v(1, :) = v_perpendicular*cos(omega*t);
v(2, :) = -v_perpendicular*sin(omega*t);
v(3, :) = v_parallel;

x(1, :) = rl*sin(omega*t) + x0(1);
x(2, :) = rl*cos(omega*t) + x0(2)-rl;
x(3, :) = v_parallel*t + x0(3);