T0 = 0; T = 20*pi;
N = 1000;
dt = (T-T0)/N;
q = 1; m = 1;
x0 = [0, 1, 0]'; v0 = [1, 0, 0]';
B0 = 1;

% Normalization
Omega = q*B0/m;
vp = sqrt(v0(1)*v0(1)+v0(2)*v0(2));
r_L = vp/Omega;
x0 = x0/r_L;
v0 = v0/vp;
T0 = T0*Omega;
T = T*Omega;
dt = dt*Omega;

B = @(x, t) B0*[0, 0, 1]';
E = @(x, t) 0.1*vp*B0*[0, 1, 0]';


[x, v, t] = larmor_motion_dimensionless_solver(E, B, x0, v0, T0, T, dt, 2);
[xt, vt, tt] = larmor_motion_analytical_2(N);

% Plot trajectory
figure;
plot(x(1, :), x(2, :), 'r:', 'LineWidth', 1.5); hold on;
plot(xt(1, :), xt(2, :), 'b-', 'LineWidth', 1.5); hold on;
legend('numerical', 'analytical');
line([0 0], ylim, 'LineWidth', 1.0, 'Color', 'k', 'HandleVisibility', 'off');
line(xlim, [0 0], 'LineWidth', 1.0, 'Color', 'k', 'HandleVisibility', 'off'); 
xlabel("Position, $x'$", 'Interpreter', 'latex');
ylabel("Position, $y'$", 'Interpreter', 'latex');
set(gca, 'FontSize', 12);
title('Trajectory');

% Plot x position
figure;
plot(t, x(1, :), 'r:', 'LineWidth', 2.0); hold on;
plot(t, xt(1, :), 'b-', 'LineWidth', 1.5); hold on;
legend('numerical', 'analytical');
line([0 0], ylim, 'LineWidth', 1.0, 'Color', 'k', 'HandleVisibility', 'off');
line(xlim, [0 0], 'LineWidth', 1.0, 'Color', 'k', 'HandleVisibility', 'off');
xlabel("Time, $t'$", 'Interpreter', 'latex');
ylabel("Position, $x'$", 'Interpreter', 'latex');
set(gca, 'FontSize', 12);
title('x-Position vs. Time');

errors = [];
% Ns = 10.^[3:6];
Ns = 10.^[3:0.5:6];
for N = Ns
    dt = (T-T0)/N;
    [x, v, t] = larmor_motion_dimensionless_solver(E, B, x0, v0, T0, T, dt, 2);
    [xt, vt, tt] = larmor_motion_analytical_2(N);
    xt = xt/r_L; vt = vt/vp; tt = tt*Omega;
    errors = [errors norm(xt(:, end)-x(:, end))];
end
figure;
loglog(Ns, errors, '*-', 'LineWidth', 2.0); grid on;
xlabel("Number of timesteps, $N$", 'Interpreter', 'latex');
ylabel("2-norm of errors in position $x$ at $t = 20\pi$", 'Interpreter', 'latex');
title("Error vs. timesteps");