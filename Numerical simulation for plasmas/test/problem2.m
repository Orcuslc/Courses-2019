T0 = 0; T = 20*pi;
N = 2000000;
dt = (T-T0)/N;
q = 1; m = 1;
x0 = [0, 1, 0]'; v0 = [1, 0, 0]';

B1 = @(x, t) [0, 0, 1]';
E2 = @(x, t) [0, 0.1, 0]';

[x, v, t] = larmor_motion_solver(E2, B1, q, m, x0, v0, T0, T, dt);
[xt, vt, tt] = larmor_motion_analytical_ExB_drift(N);

% Plot trajectory
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

% Plot x position
figure;
plot(t, x(1, :), 'r:', 'LineWidth', 2.0); hold on;
plot(t, xt(1, :), 'b-', 'LineWidth', 1.5); hold on;
legend('numerical', 'analytical');
line([0 0], ylim, 'LineWidth', 1.0, 'Color', 'k', 'HandleVisibility', 'off');
line(xlim, [0 0], 'LineWidth', 1.0, 'Color', 'k', 'HandleVisibility', 'off');
xlabel('Time, $t$', 'Interpreter', 'latex');
ylabel('Position, $x$', 'Interpreter', 'latex');
set(gca, 'FontSize', 12);
title('x-Position vs. Time');

% errors = [];
% Ns = 10.^[3:6];
% % Ns = 10.^[3:0.5:6];
% for N = Ns
%     dt = (T-T0)/N;
%     [x, v, t] = larmor_motion_solver(E2, B1, q, m, x0, v0, T0, T, dt, 1);
%     [xt, vt, tt] = larmor_motion_analytical_2(N);
%     errors = [errors norm(xt(:, end)-x(:, end))];
% end
% figure;
% loglog(Ns, errors, '*-', 'LineWidth', 2.0); grid on;
% xlabel("Number of timesteps, $N$", 'Interpreter', 'latex');
% ylabel("2-norm of errors in position $x$ at $t = 20\pi$", 'Interpreter', 'latex');
% title("Error vs. timesteps");