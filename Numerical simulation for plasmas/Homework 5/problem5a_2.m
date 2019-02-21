
vpars = [];
vpreps = [];

taos_nu = [];
% collision frequency
[t1, vpars, vpreps] = monte_carlo_collisions(0.001, vpars, vpreps);
[t2, vpars, vpreps] = monte_carlo_collisions(0.01, vpars, vpreps);
[t3, vpars, vpreps] = monte_carlo_collisions(0.1, vpars, vpreps);
[t4, vpars, vpreps] = monte_carlo_collisions(1, vpars, vpreps);

figure;
semilogx([0.001, 0.01, 0.1, 1], [mean(t1), mean(t2), mean(t3), mean(t4)], '-*'); hold on;
errorbar([0.001, 0.01, 0.1, 1], [mean(t1), mean(t2), mean(t3), mean(t4)], [std(t1), std(t2), std(t3), std(t4)]);
xlim([1e-4 10]);
xlabel("$\nu/\Omega$", "Interpreter", "latex");
ylabel("$\tau$", "Interpreter", "latex");

figure;
plot(vpars, vpreps, '.');
xlabel("$v_\parallel$", "Interpreter", "latex"); 
ylabel("$v_\perp$", "Interpreter", "latex");


function [taos, vpars, vpreps] = monte_carlo_collisions(nu, vpars, vpreps)
% taos: confinement time;
    taos = [];
    
    for m = 1:50
        m
        % dimensionless arguments
        B = @B_4a;
        E = @(x, t) [0, 0, 0]';

        x0 = [0, 0, 0]';
        v0 = [0, 1, 1]';
        T_end = 100*pi;
        dT = T_end/5000;
        Ts = 0:dT:T_end;
        L = 10;
        
        % confinement time
        tao = 0;
        
        for i = 1:5000
            [x, v, t] = larmor_motion_dimensionless_solver(E, B, x0, v0, Ts(i), Ts(i+1), 1e-3, 'RK45', 1e-3);

            % check if particle escaped the magnetic mirror
            index = find(abs(x(3, :)) > L/2);
            if(~isempty(index)) % escaped
                tao = t(index(1));
                taos = [taos, tao];
                [vpar, vprep] = vproj(x(:, index(1)), v(:, index(1)), t(index(1)));
                vpars = [vpars, vpar];
                vpreps = [vpreps, vprep];
                break;
            end

            % not escaped
            x0 = x(:, end); 
            v0 = v(:, end);

            % check if collision occured
            p_collision = 1-exp(-nu*dT);
            c = rand();
            if(p_collision > c) % occured
                v0 = small_angle_collision(v0, pi/18);
            end
        end
    end
end

function [vpar, vprep] = vproj(x, v, t)
B = B_4a(x, t);
vpar = v'*B/norm(B)^2*B;
vprep = v - vpar;
vpar = norm(vpar);
vprep = norm(vprep);
end
