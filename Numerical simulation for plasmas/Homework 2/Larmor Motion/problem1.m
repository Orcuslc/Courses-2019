T0 = 0; T = 20*pi;
dt = 1e-3;
% dt = 4e-2;
q = 1; m = 1;
x0 = [0, 1, 0]'; v0 = [1, 0, 0]';

[x, v, t] = larmor_motion_solver(@E1, @B1, q, m, x0, v0, T0, T, dt);

% Plot trajectory
figure;
plot(x(1, 1), x(2, 1), '*', 'LineWidth', 6.0); hold on;
plot(x(1, :), x(2, :), 'b-', 'LineWidth', 1.5); hold on;
plot(x(1, :), x(2, :), 'r:', 'LineWidth', 1.5); hold on;
line([0 0], ylim, 'LineWidth', 1.0, 'Color', 'k');
line(xlim, [0 0], 'LineWidth', 1.0, 'Color', 'k'); 
xlabel('Position, $x$', 'Interpreter', 'latex');
ylabel('Position, $y$', 'Interpreter', 'latex');
set(gca, 'FontSize', 12);
title('Trajectory');

% Plot x position
figure;
plot(t, x(1, :), 'b-', 'LineWidth', 2.0); hold on;
plot(t, x(1, :), 'r:', 'LineWidth', 2.0); hold on;
line([0 0], ylim, 'LineWidth', 1.0, 'Color', 'k');
line(xlim, [0 0], 'LineWidth', 1.0, 'Color', 'k');
xlabel('Time, $t$', 'Interpreter', 'latex');
ylabel('Position, $x$', 'Interpreter', 'latex');
set(gca, 'FontSize', 12);
title('x-Position vs. Time');
