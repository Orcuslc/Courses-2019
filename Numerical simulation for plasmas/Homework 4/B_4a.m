function B = B_4a(x, t)
% Normalized magnetic field
% x, t: normalized parameter

B = zeros(3, 1);

% Normalized parameters
L = 10;
B00 = 10;
Rm = 4;
delta_Bz = B00*(Rm-1);

z = x(3); y = x(2); x = x(1);

B(1) = -pi*x/(2*L)*delta_Bz*sin(2*pi*z/L);
B(2) = -pi*y/(2*L)*delta_Bz*sin(2*pi*z/L);
B(3) = B00 + delta_Bz/2*(1-cos(2*pi*z/L));