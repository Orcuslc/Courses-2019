function [x, v, t] = larmor_motion_dimensionless_solver(E, B, x0, v0, T0, T, dt, order, rel_tol)
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
% order: 1 (Euler)
%     or 2 (Leapfrog)
%     or AB2 (Adams-Bashforth, second order)
%     or AB3 (Adams-Bashforth, third order)
%     or RK45 (Runge-Kutta 4th order, with 5th order correction)
% rel_tol: for RK45. default: 1e-6;

% Output: 
% x, v: 3xL matrix, each column is a solution
% t: 1xL array, the timestamps

t = T0:dt:T;
x = zeros(3, length(t));
v = zeros(3, length(t));
x(:, 1) = x0; v(:, 1) = v0;

drift = @(v, x, t) E(x, t)+cross(v, B(x, t));

if nargin == 7
    order = 1; % default: Euler
    rel_tol = 1e-6;
elseif nargin == 8
    rel_tol = 1e-6;
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
    case 'AB2'
        % Adams-Bashforth, second order
        % the first step is to use Euler method, then use AB2
        x(:, 2) = x(:, 1) + v(:, 1)*dt;
        v(:, 2) = v(:, 1) + drift(v(:, 1), x(:, 1), t(1))*dt;
        for i = 2:(length(t)-1)
            x(:, i+1) = x(:, i) + 0.5*(3*v(:, i) - v(:, i-1))*dt;
            v(:, i+1) = v(:, i) + 0.5*(3*drift(v(:, i), x(:, i), t(i)) - drift(v(:, i-1), x(:, i-1), t(i-1)))*dt;
        end
    case 'AB3'
        % Adams-Bashforth, third order
        % use Euler for x0->x0.25, then use leapfrog to get x0.5, x1, x2.
        
        % x0 - >x0.25
        x1 = x0 + v0*0.25*dt;
        v1 = x0 + drift(v0, x0, T0)*0.25*dt;
        
        % x0, x0.25 -> x0.5
        x2 = x0 + v1*0.5*dt;
        v2 = v0 + drift(v1, x1, T0+0.25*dt)*0.5*dt;
        
        % x0, x0.5 -> x1
        x(:, 2) = x0 + v2*dt;
        v(:, 2) = v0 + drift(v2, x2, T0+0.5*dt)*dt;
        
        % x0, x1 -> x2
        x(:, 3) = x0 + v(:, 2)*2*dt;
        v(:, 3) = v0 + drift(v(:, 2), x(:, 2), t(2))*2*dt;  
        
        for i = 3:(length(t)-1)
            x(:, i+1) = x(:, i) + (23*v(:, i)-16*v(:, i-1)+5*v(:, i-2))*dt/12;
            v(:, i+1) = v(:, i) + (23*drift(v(:, i), x(:, i), t(i))-16*drift(v(:, i-1), x(:, i-1), t(i-1))+5*drift(v(:, i-2), x(:, i-2), t(i-2)))*dt/12;
        end
    case 'RK45'
        % Runge-Kutta 4th order, with 5th order correction
        % use `ode45` in matlab for integration
        % the equation becomes
        % d^2x/dt^2 = (E+dx/dt \cross B);
        
        y0 = [x0', v0'];
        func = @(t, x) larmor_motion_ode_func(B, E, t, x);
        options = odeset('RelTol', rel_tol);
        [t, y] = ode45(func, [T0 T], y0, options);
        x = y(:, 1:3)';
        v = y(:, 4:6)';
        t = t';
end
        