function [xt, vt, t] = larmor_motion_analytical_1(N)
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
E0=[0. 0. 0.]';   % Initial Electric Field (Column vector)
Bmag=norm(B0);    % Magnitude of B0 (2-norm)

% Compute timestep
dt=(tspan(end)-tspan(1))/N;
%Calculate row vector of times
t=(tspan(1):dt:tspan(end))';
%t=t';
%Compute total number of times
nt=size(t,1);

%========================================================================
%Compute Analytic Solution
%========================================================================
xt=zeros(3,nt);
vt=zeros(3,nt);
rl=m*sqrt(v0(1)*v0(1)+v0(2)*v0(2))/(q*B0(3));
om=q*B0(3)/m;
xt(1,:)=rl*sin(om*t) + x0(1);
xt(2,:)=rl*cos(om*t) + (x0(2)-rl);