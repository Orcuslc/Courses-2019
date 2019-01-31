function B = B_3a_2(x, t)
% DIPOLE MAGNETIC FIELD
% Dipole moment (sets strength of magnetic field)
% NOTE: This includes [ mu_0/(4 pi B_0) ] factor (dimensionless norm)
dm=-100.0;  

B = zeros(3, 1);

r=sqrt(x(1).*x(1)+x(2).*x(2)+x(3).*x(3)); % Spherical radius
cr=sqrt(x(1).*x(1)+x(2).*x(2)); % Cylindrical radius
B(1)=3.0*x(1).*x(3)./r.^2;
B(2)=3.0*x(2).*x(3)./r.^2;
B(3)=2.0*x(3).*x(3)./r.^2 - cr.*cr./r.^2;

%Add 1/r^3 dependence and dipole moment
B=dm*B./r.^3;