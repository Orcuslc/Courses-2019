% EXAMPLE Matlab Code for Analytical Solution of Single Particle Motion
% Input x0, v0, B0, E0,
% Compute x and v as a function of time
%-----------------------------------------------------------------------
% Professor Gregory Howes
% PHYS:5905 Numerical Simulation of Plasmas
% 16 JAN 2019
'======================================================================='

%Setting TimeSpan and Initial Conditions
tspan= [0. 10.*2.*pi]';
nn=1000; %Total number of timesteps to advance
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
% E0=[0. 0. 0.]';   % Initial Electric Field (Column vector)
E0 = [0, 0.1, 0]';
Bmag=norm(B0);    % Magnitude of B0 (2-norm)

% Compute timestep
dt=(tspan(end)-tspan(1))/nn;
%Calculate row vector of times
t=(tspan(1):dt:tspan(end))';
%t=t';
%Compute total number of times
nt=size(t,1)

%========================================================================
%Compute Analytic Solution
%========================================================================
xt=zeros(3,nt);
vt=zeros(3,nt);
rl=m*sqrt(v0(1)*v0(1)+v0(2)*v0(2))/(q*B0(3));
om=q*B0(3)/m;
xt(1,:)=rl*sin(om*t) + x0(1);
xt(2,:)=rl*cos(om*t) + (x0(2)-rl);

%Note that I could also compute Analytical solution using a Loop
xt2=zeros(3,nt);
vt2=zeros(3,nt);
% Loop over timesteps
for i= 1:nt
xt2(1,i)=rl*sin(om*t(i)) + x0(1);
xt2(2,i)=rl*cos(om*t(i)) + (x0(2)-rl);
end

%========================================================================
% PLOT FIGURES
%========================================================================
%Get Screen size
scrsz = get(0,'ScreenSize');

%Set Figure Size
h1=figure('Position',[1 scrsz(4) 1*scrsz(3)/2 2*scrsz(4)/3]);
% Plot (x,y) Trajectory
plot(xt(1,:),xt(2,:),'LineWidth',6.0,'Color','b','LineStyle','-')
hold on
plot(xt2(1,:),xt2(2,:),'LineWidth',2.0,'Color','r','LineStyle',':');
% Plot Axes
xL = xlim;
yL = ylim;
line([0 0], yL,'LineWidth',1.0,'Color','k');  %y-axis
line(xL, [0 0],'LineWidth',1.0,'Color','k');  %x-axis
%Plot Labels
xlabel('Position, $x$','Interpreter','latex','FontName','TimesNewRoman','FontSize',20,'FontWeight','bold')
ylabel('Position, $y$','Interpreter','latex','FontName','TimesNewRoman','FontSize',20,'FontWeight','bold')
set(gca,'FontSize',16,'FontName','TimesNewRoman','FontWeight','bold','LineWidth',2)
title('Trajectory');


%Set Figure Size
% Plot x(t) vs. t
h2=figure('Position',[1*scrsz(3)/2 scrsz(4) 1*scrsz(3)/2 2*scrsz(4)/3]);
plot(t(:),xt(1,:),'LineWidth',2.0,'Color','b','LineStyle','-');
hold on
plot(t(:),xt2(1,:),'LineWidth',2.0,'Color','r','LineStyle',':');
hold off
% Plot Axes
xL = xlim;
yL = ylim;
line([0 0], yL,'LineWidth',1.0,'Color','k');  %y-axis
line(xL, [0 0],'LineWidth',1.0,'Color','k');  %x-axis
%Plot Labels
xlabel('Time, $t$','Interpreter','latex','FontName','TimesNewRoman','FontSize',20,'FontWeight','bold')
ylabel('Position, $x$','Interpreter','latex','FontName','TimesNewRoman','FontSize',20,'FontWeight','bold')
set(gca,'FontSize',16,'FontName','TimesNewRoman','FontWeight','bold','LineWidth',2)
title('x-Position vs. Time');


