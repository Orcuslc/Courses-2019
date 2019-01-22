function [xt, vt, t] = larmor_motion_analytical_dimensionless(N)
% N: number of timesteps

%Setting TimeSpan and Initial Conditions
tspan= [0. 10.*2.*pi]';
%Initial Conditions
x0=zeros(3,1);
v0=zeros(3,1);
x0 = [0. 1. 0. ]';  % Initial position (Column vector)
v0 = [1. 0. 0. ]';  % Initial velocity (Column vector)

%Parameters
q=1.;  %Charge
m=1.;  %Mass

%Electromagnetic Field Input
B0=[0. 0. 1.]';   % Initial Magnetic Field (Column vector)
E0=[0. 0.1 0.]';   % Initial Electric Field (Column vector)
Bmag=norm(B0);    % Magnitude of B0 (2-norm)

% Compute timestep
dt=(tspan(end)-tspan(1))/N;
%Calculate row vector of times
t=(tspan(1):dt:tspan(end))';
%t=t';
%Compute total number of times
nt=size(t,1);

% Normalization
vp = sqrt(v0(1)^2+v0(2)^2);
v0 = v0/vp;
E0 = E/(vp*Bmag);
B0 = B/Bmag;
Omega = q*Bmag/m;
t = t*Omega;
r_L = vp/Omega;
x0 = x0/r_L;

%========================================================================
%Compute Analytic Solution
%========================================================================
xt=zeros(3,nt);
vt=zeros(3,nt);

om=q*B0(3)/m;
vp = sqrt(v0(1)*v0(1)+v0(2)*v0(2)); % perpendicular part of velocity
drift = cross(E0, B0)/sum(B0.^2);

vt(1, :) = vp*cos(om*t);
vt(2, :) = -vp*sin(om*t);
vt(3, :) = (q/m*E0(3)*t+v0(3));
vt = vt + drift;

xt(1, :) = x0(1) + vp/om*sin(om*t) + drift(1)*t;
xt(2, :) = x0(2)-vp/om + vp/om*cos(om*t) + drift(2)*t;
xt(3, :) = x0(3) + 1/2*q/m*E0(3)*(t.^2) + v0(3)*t;

% rl=m*sqrt(v0(1)*v0(1)+v0(2)*v0(2))/(q*B0(3));
% xt(1,:)=rl*sin(om*t) + x0(1);
% xt(2,:)=rl*cos(om*t) + (x0(2)-rl);