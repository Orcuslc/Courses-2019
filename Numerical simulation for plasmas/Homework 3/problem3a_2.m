T0 = 0; 
% T = 30*pi;
T = 60*pi;
N = 10000;
dt = (T-T0)/N;

% Normalized initial conditions
x0 = [5., 0., 0.]';
v0 = [0., -1., 0.]';


B = @B_3a_2;
E = @(x, t) [0, 0, 0]';


[x, v, t] = larmor_motion_dimensionless_solver(E, B, x0, v0, T0, T, dt, 'AB3');

% Plot trajectory
figure;
plot(x(1, :), x(2, :), 'r:', 'LineWidth', 1.5); hold on;
legend('numerical');
line([0 0], ylim, 'LineWidth', 1.0, 'Color', 'k', 'HandleVisibility', 'off');
line(xlim, [0 0], 'LineWidth', 1.0, 'Color', 'k', 'HandleVisibility', 'off'); 
xlabel("Position, $x'$", 'Interpreter', 'latex');
ylabel("Position, $y'$", 'Interpreter', 'latex');
set(gca, 'FontSize', 12);
title('Trajectory');