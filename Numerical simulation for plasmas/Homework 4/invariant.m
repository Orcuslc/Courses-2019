function [magnetic_moment, magnetic_magnitude, parallel_kinetic_energy, kinetic_energy, perpendicular_velocity, parallel_velocity] ...
    = invariant(x, v, t)

% notice: all the energy are of unit mass, i.e., m = 1

B = B_4a(x, t);
magnetic_magnitude = norm(B);

% angle between velocity and magnetic field
alpha = acos(v'*B/(norm(v)*norm(B)));

% components of velocity
perpendicular_velocity = v*sin(alpha);
parallel_velocity = v*cos(alpha);

magnetic_moment = norm(perpendicular_velocity)^2/(2*magnetic_magnitude);
parallel_kinetic_energy = 1/2*norm(parallel_velocity)^2;
kinetic_energy = 1/2*norm(v)^2;