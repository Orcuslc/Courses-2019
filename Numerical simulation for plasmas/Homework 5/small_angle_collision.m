function y0 = small_angle_collision(y0, sigma)
% Compute a small angle collision
% y0: [x, y, z, vx, vy, vz]
% sigma: standard deviation of rotation angle theta

% y: after rotation

[phi1, elev1, radius] = cart2sph(y0(4), y0(5), y0(6));
theta1 = elev1+pi/2;

% rotate to z-axis
R_z1 = [cos(-phi1), sin(-phi1), 0;
    -sin(-phi1), cos(-phi1), 0;
    0, 0, 1];
R_y1 = [cos(-theta1), 0, sin(-theta1);
    0, 1, 0;
    -sin(-theta1), 0, cos(-theta1)];

% sigma: the sd of theta
theta = sigma*randn();

% phi: unif([0, 2*pi]);
phi = 2*pi*rand();

% rotate matrix w.r.t. sigma and phi
R1 = [cos(theta), 0, -sin(theta);
      0, 1, 0;
      sin(theta), 0, cos(theta)];
R2 = [cos(phi), sin(phi), 0;
      -sin(phi), cos(phi), 0;
      0, 0, 1];

% rotate back
R3 = [cos(-theta1), 0, sin(-theta1);
      0, 1, 0;
      -sin(-theta1), 0, cos(theta1)];
R4 = [cos(phi1), sin(phi1), 0;
      -sin(phi1), cos(phi1), 0;
      0, 0, 1];
 
y0(4:6) = R4*R3*R2*R1*R_y1*R_z1*y0(4:6)';
