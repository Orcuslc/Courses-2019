function dy = larmor_motion_ode_func(B, E, t, y)
% second-order ode function for larmor motion d^2x/dt^2 = (E+dx/dt \cross B);
% y: row vector, [x, y, z, vx, vy, vz];
% dy: row vector, [dx/dt, dy/dt, dz/dt, dvx/dt, dvy/dt, dvz/dt];

dy = zeros(size(y));
dy(1:3) = y(4:6);
dy(4:6) = (E(y(1:3)', t) + cross(y(4:6)', B(y(1:3)', t))')';