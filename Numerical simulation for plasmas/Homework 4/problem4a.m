T0 = 0; 
T = 10*pi;
N = 20000;
dt = (T-T0)/N;

x0 = [0.1, 0, 8]'; 
v0 = [0, 1, 0]';

B = @B_4a;
E = @(x, t) [0, 0, 0]';

[x, v, t] = larmor_motion_dimensionless_solver(E, B, x0, v0, T0, T, dt, 'AB3');

% Plot trajectory
figure;
plot(x(3, :), x(1, :), 'r:', 'LineWidth', 1.5); hold on;
xlabel("Position, $z'$", 'Interpreter', 'latex');
ylabel("Position, $x'$", 'Interpreter', 'latex');
set(gca, 'FontSize', 12);
title('Trajectory');

% 3D trajectory
figure;
plot3(x(1, :), x(2, :), x(3, :), 'r:');
xlabel("Position, $x'$", 'Interpreter', 'latex');
ylabel("Position, $y'$", 'Interpreter', 'latex');
zlabel("Position, $z'$", 'Interpreter', 'latex');
axis([-0.5, 0.5, -0.5, 0.5, -10, 10]);
title("3D trajectory");

magnetic_moments = zeros(size(t));
magnetic_magnitudes = zeros(size(t));
parallel_kinetic_energies = zeros(size(t));
total_kinetic_energies = zeros(size(t));
perpendicular_velocities = zeros(size(t));

for i = 1:size(t, 2)
    [magnetic_moment, magnetic_magnitude, parallel_kinetic_energy, kinetic_energy, perpendicular_velocity, parallel_velocity] ...
        = invariant(x(:, i), v(:, i), t(i));
    magnetic_moments(i) = magnetic_moment;
    magnetic_magnitudes(i) = magnetic_magnitude;
    parallel_kinetic_energies(i) = parallel_kinetic_energy;
    total_kinetic_energies(i) = kinetic_energy;
    perpendicular_velocities(i) = norm(perpendicular_velocity);
end

% magnetic moment and perpendicular velocity
figure;
plot(t, magnetic_moments/magnetic_moments(1), ':', 'LineWidth', 2.5); hold on;
plot(t, perpendicular_velocities.^2/perpendicular_velocities(1)^2, '--', 'LineWidth', 3.5); hold on;
plot(t, magnetic_magnitudes/magnetic_magnitudes(1), '-', 'LineWidth', 1.5); hold on;
legend("$\mu$", "$v_\perp^2$", "$B$", "Interpreter", "latex");
xlabel("$t'$", "Interpreter", "latex");
ylabel("percentage");

% kinetic energy
figure;
plot(t, total_kinetic_energies/total_kinetic_energies(1), ':', 'LineWidth', 2.5); hold on;
plot(t, parallel_kinetic_energies/total_kinetic_energies(1), '--', 'LineWidth', 2.5); hold on;
plot(t, magnetic_moments.*magnetic_magnitudes/total_kinetic_energies(1), '-', 'LineWidth', 1.5); hold on;
legend("$E$", "$E_\parallel$", "$\mu B$", "Interpreter", "latex");
xlabel("$t'$", "Interpreter", "latex");
ylabel("percentage");