function B = magfield(t,x)
% MAGNETIC Calculates magnetic field B(t,x)
% Input values are scalar time t and vector position x(3)= (x,y,z)
% Output values are column vector B(3)=(Bx,By,Bz)
%
%  FIELD Geometry choice (geom):
%    1 Straight Field
%    2 Dipole Field
  geom=2;

% Initialize column vector
  B=zeros(3,1);

switch(geom)
 case 1
% Straight Magnetic Field
B(1)=0. ;
B(2)=0. ;
B(3)=1.0 ;
case 2
%DIPOLE MAGNETIC FIELD
%Dipole moment (sets strength of magnetic field)
%    NOTE: This includes [ mu_0/(4 pi B_0) ] factor (dimensionless norm)
dm=-1000.0;
  
r=sqrt(x(1).*x(1)+x(2).*x(2)+x(3).*x(3)); % Spherical radius
cr=sqrt(x(1).*x(1)+x(2).*x(2)); % Cylindrical radius
B(1)=3.0*x(1).*x(3)./r.^2;
B(2)=3.0*x(2).*x(3)./r.^2;
B(3)=2.0*x(3).*x(3)./r.^2 - cr.*cr./r.^2;

%Add 1/r^3 dependence and dipole moment
B=dm*B./r.^3;

end
