% Problem 1

figure;
for i = 1:1000
    y0 = [1, 0, 0]'; % [x, v]
    y1 = small_angle_collision(y0, pi/18);
    plot3(y1(1), y1(2), y1(3), '.'); hold on;
end
axis([-1, 1, -1, 1, -1, 1]);
xlabel('x'); ylabel('y'); zlabel('z');
title('$\sigma(\theta) = \pi/18$', 'Interpreter', 'latex');
grid on;

figure;
for i = 1:1000
    y0 = [1, 0, 0]'; % [x, v]
    y1 = small_angle_collision(y0, pi/180);
    plot3(y1(1), y1(2), y1(3), '.'); hold on;
end
axis([-1, 1, -1, 1, -1, 1]);
xlabel('x'); ylabel('y'); zlabel('z');
title('$\sigma(\theta) = \pi/180$', 'Interpreter', 'latex');
grid on;