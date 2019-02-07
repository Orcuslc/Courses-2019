T0 = 0; T = 20*pi;
N = 10000;
dt = (T-T0)/N;

x0 = [0, 1, 0]'; v0 = [1, 0, 0]';

B = @(x, t) [0, 0, 1]';
E = @(x, t) [0, 0.1, 0]';

[x1, v1, t1] = larmor_motion_dimensionless_solver(E, B, x0, v0, T0, T, dt, 'AB3');
[xt, vt, tt] = larmor_motion_analytical_ExB_drift(N);

error_AB3 = norm(x1(:, end)-xt(:, end))
[x2, v2, t2] = larmor_motion_dimensionless_solver(E, B, x0, v0, T0, T, dt, 'RK45', error_AB3);
size(t2, 2)

[x, v, t] = larmor_motion_dimensionless_solver(E, B, x0, v0, T0, T, dt, 'RK45', 4e-3);
size(t, 2)
error_RK45 = norm(x(:, end)-xt(:, end))/norm(xt(:, end))

% % Plot trajectory
figure;
plot(x(1, :), x(2, :), 'r:', 'LineWidth', 1.5); hold on;
plot(xt(1, :), xt(2, :), 'b-', 'LineWidth', 1.5); hold on;
legend('numerical', 'analytical');
line([0 0], ylim, 'LineWidth', 1.0, 'Color', 'k', 'HandleVisibility', 'off');
line(xlim, [0 0], 'LineWidth', 1.0, 'Color', 'k', 'HandleVisibility', 'off'); 
xlabel('Position, $x$', 'Interpreter', 'latex');
ylabel('Position, $y$', 'Interpreter', 'latex');
set(gca, 'FontSize', 12);
title('Trajectory');


[x2, v2, t2] = larmor_motion_dimensionless_solver(E, B, x0, v0, T0, T, dt, 'RK45', 1e-6);
size(t2, 2)
error_RK45 = norm(x2(:, end)-xt(:, end))/norm(xt(:, end))


% Energy
N = 2208;
dt = (T-T0)/N;
[x, v, t] = larmor_motion_dimensionless_solver(E, B, x0, v0, T0, T, dt, 'AB3');

total_kinetic_energies = zeros(size(t));
for i = 1:size(t, 2)
    [magnetic_moment, magnetic_magnitude, parallel_kinetic_energy, kinetic_energy, perpendicular_velocity, parallel_velocity] ...
        = invariant(x(:, i), v(:, i), t(i));
    total_kinetic_energies(i) = kinetic_energy;
end

1-total_kinetic_energies(end)/total_kinetic_energies(1)


[x, v, t] = larmor_motion_dimensionless_solver(E, B, x0, v0, T0, T, dt, 'RK45', 1e-3);
total_kinetic_energies = zeros(size(t));
for i = 1:size(t, 2)
    [magnetic_moment, magnetic_magnitude, parallel_kinetic_energy, kinetic_energy, perpendicular_velocity, parallel_velocity] ...
        = invariant(x(:, i), v(:, i), t(i));
    total_kinetic_energies(i) = kinetic_energy;
end
size(t, 2)
1-total_kinetic_energies(end)/total_kinetic_energies(1)
