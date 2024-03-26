% Higher Order ODE: Example with a second-order ODE
% y'' + 2y' + y = 0
% Convert to a system of first-order ODEs

% Define the ODE system
ode_fun = @(t, Y) [Y(2); -2*Y(2) - Y(1)];

% Time span
tspan = [0, 10];

% Initial conditions for y and y'
initial_conditions = [1, 0];

% Solve the system
[t, Y] = ode45(ode_fun, tspan, initial_conditions);

% Plot the results
figure;
plot(t, Y(:, 1), '-o', 'DisplayName', 'y(t)');
hold on;
plot(t, Y(:, 2), '-o', 'DisplayName', "y'(t)");
xlabel('Time (t)');
ylabel('Solution');
title('Higher Order ODE: y'''' + 2y'' + y = 0');
legend;
grid on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Pendulum Simulation

% Define the ODE system for a simple pendulum
pendulum_ode = @(t, Y) [Y(2); -sin(Y(1))];

% Time span
tspan_pendulum = [0, 10];

% Initial conditions for angle and angular velocity
initial_conditions_pendulum = [pi/4, 0];

% Solve the system
[t_pendulum, Y_pendulum] = ode45(pendulum_ode, tspan_pendulum, initial_conditions_pendulum);

% Convert polar coordinates to Cartesian for plotting
x_pendulum = sin(Y_pendulum(:, 1));
y_pendulum = -cos(Y_pendulum(:, 1));

% Plot the results
figure;
plot(t_pendulum, x_pendulum, '-o', 'DisplayName', 'x(t)');
hold on;
plot(t_pendulum, y_pendulum, '-o', 'DisplayName', 'y(t)');
xlabel('Time (t)');
ylabel('Position');
title('Pendulum Simulation');
legend;
grid on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Orbital Mechanics Simulation

% Define the ODE system for orbital mechanics
orbital_ode = @(t, Y) [-Y(1) / norm(Y(1:2))^3; -Y(2) / norm(Y(1:2))^3; Y(3); Y(4)];

% Time span
tspan_orbital = [0, 100];

% Initial conditions for position and velocity
initial_conditions_orbital = [1, 0, 0, 1];

% Solve the system
[t_orbital, Y_orbital] = ode45(orbital_ode, tspan_orbital, initial_conditions_orbital);

% Plot the results
figure;
plot(Y_orbital(:, 1), Y_orbital(:, 2), '-o', 'DisplayName', 'Orbit');
xlabel('X-axis');
ylabel('Y-axis');
title('Orbital Mechanics Simulation');
legend;
grid on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Program for Example 6.13: Euler's Method for a System of Equations
% Input: interval inter, initial vector y0, number of steps n
% Output: time steps t, solution y

function [t, y] = euler_system_example(inter, y0, n)
    t(1) = inter(1);
    y(1, :) = y0;
    h = (inter(2) - inter(1)) / n;

    for i = 1:n
        t(i+1) = t(i) + h;
        y(i+1, :) = eulerstep_system(t(i), y(i, :), h);
    end

    % Plot the results
    figure;
    plot(t, y(:, 1), '-o', 'DisplayName', 'y1 - Euler Method');
    hold on;
    plot(t, y(:, 2), '-s', 'DisplayName', 'y2 - Euler Method');
    
    % Include the correct solution for comparison
    exact_solution = @(t) [t .* exp(-2 * t); exp(-t)];
    exact_values = exact_solution(t);
    plot(t, exact_values(1, :), '--', 'DisplayName', 'y1 - Exact Solution');
    plot(t, exact_values(2, :), '--', 'DisplayName', 'y2 - Exact Solution');

    xlabel('t');
    ylabel('Approximate Solution / Exact Solution');
    title('Euler Method for System of Equations (Example 6.13)');
    legend('show');
end

function y = eulerstep_system(t, y, h)
    % One step of the Euler Method for a system
    % Input: current time t, current vector y, step size h
    % Output: the approximate solution vector at time t + h
    y = y + h * ydot_system(t, y)';
end

function z = ydot_system(t, y)
    % Right-hand side of the system of differential equations
    z(1) = y(2)^2 - 2 * y(1);
    z(2) = y(1) - y(2) - t * y(2)^2;
end


euler_system_example([0 1], [0 1], 10);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Program for Example 6.13: Simulation of a Pendulum using the Trapezoid Method
% Input: time interval inter, initial values ic = [y1(0) y2(0)], number of steps n
% Output: Animated visualization of the pendulum motion

function pendulum_simulation(inter, ic, n)
    h = (inter(2) - inter(1)) / n; % time step
    y(1, :) = ic; % initial conditions
    t(1) = inter(1);

    % Set up the plot
    figure;
    set(gca, 'xlim', [-1.2 1.2], 'ylim', [-1.2 1.2], ...
        'XTick', [-1 0 1], 'YTick', [-1 0 1]);
    title('Simulation of Pendulum Motion');
    xlabel('x');
    ylabel('y');
    hold on;

    bob = animatedline('color', 'r', 'Marker', '.', 'markersize', 40);
    rod = animatedline('color', 'b', 'LineStyle', '-', 'LineWidth', 3);
    axis square;

    for k = 1:n
        t(k+1) = t(k) + h;
        y(k+1, :) = trapstep(t(k), y(k, :), h);

        % Update pendulum position
        x_bob = sin(y(k+1, 1));
        y_bob = -cos(y(k+1, 1));
        x_rod = [0 x_bob];
        y_rod = [0 y_bob];

        % Update animated points
        clearpoints(bob);
        addpoints(bob, x_bob, y_bob);
        clearpoints(rod);
        addpoints(rod, x_rod, y_rod);

        drawnow; % Update the plot
        pause(h);
    end
end

function y = trapstep(t, x, h)
    % One step of the Trapezoid Method for the pendulum equations
    z1 = ydot(t, x);
    g = x + h * z1;
    z2 = ydot(t + h, g);
    y = x + h * (z1 + z2) / 2;
end

function z = ydot(t, y)
    % Right-hand side of the pendulum equations
    g = 9.81;
    length = 1;
    z(1) = y(2);
    z(2) = -(g/length) * sin(y(1));
end

pendulum_simulation([0 10], [pi/2 0], 200);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Program for Example 6.13: Simulation of an Orbiting Satellite using the Trapezoid Method
% Input: time interval inter, initial conditions
% ic = [x0 vx0 y0 vy0], x position, x velocity, y position, y velocity,
% number of steps n, steps per point plotted p
% Calls a one-step method such as trapstep.m
% Example usage: orbit([0 100],[0 1 2 0],10000,5)

function z = orbit(inter, ic, n, p)
    h = (inter(2) - inter(1)) / n; % time step
    x0 = ic(1); vx0 = ic(2); y0 = ic(3); vy0 = ic(4); % grab initial conditions
    y(1, :) = [x0 vx0 y0 vy0]; t(1) = inter(1); % build y vector

    % Set up the plot
    figure;
    set(gca, 'XLim', [-5 5], 'YLim', [-5 5], ...
        'XTick', [-5 0 5], 'YTick', [-5 0 5]);
    sun = animatedline('color', 'y', 'Marker', '.', 'markersize', 50);
    addpoints(sun, 0, 0);
    head = animatedline('color', 'r', 'Marker', '.', 'markersize', 35);
    tail = animatedline('color', 'b', 'LineStyle', '-');
    
    for k = 1:n/p
        for i = 1:p
            t(i+1) = t(i) + h;
            y(i+1, :) = trapstep(t(i), y(i, :), h);
        end
        y(1, :) = y(p+1, :);
        t(1) = t(p+1);

        % Update plot points
        clearpoints(head);
        addpoints(head, y(1, 1), y(1, 3));
        addpoints(tail, y(1, 1), y(1, 3));
        drawnow;
    end
end

function y = trapstep(t, x, h)
    % One step of the Trapezoid Method for the orbital equations
    z1 = ydot(t, x);
    g = x + h * z1;
    z2 = ydot(t + h, g);
    y = x + h * (z1 + z2) / 2;
end

function z = ydot(t, x)
    % Right-hand side of the orbital equations
    m2 = 3; g = 1; mg2 = m2 * g; px2 = 0; py2 = 0;
    
    px1 = x(1); py1 = x(3); vx1 = x(2); vy1 = x(4);
    dist = sqrt((px2 - px1)^2 + (py2 - py1)^2);
    z = zeros(1, 4);
    
    z(1) = vx1;
    z(2) = (mg2 * (px2 - px1)) / (dist^3);
    z(3) = vy1;
    z(4) = (mg2 * (py2 - py1)) / (dist^3);
end

orbit([0 100], [0 1 2 0], 10000, 5);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Program for Example 6.18: Runge-Kutta of order four for ODE
% Example usage: rungeKuttaExample([0 1], 1, 5);

function rungeKuttaExample(inter, y0, steps)
    fprintf('steps\tstep size h\terror at t = 1\n');
    
    for n = steps
        h = (inter(2) - inter(1)) / n;
        y_approx = rungeKutta(inter, y0, n);
        exact_solution = analyticalSolution(inter(2));
        error_at_t1 = abs(y_approx - exact_solution(inter(2)));
        
        fprintf('%d\t%.5f\t%.4e\n', n, h, error_at_t1);
    end
end

function y = rungeKutta(inter, y0, n)
    h = (inter(2) - inter(1)) / n;
    t = inter(1):h:inter(2);
    y = zeros(1, n+1);
    y(1) = y0;

    for i = 1:n
        k1 = h * f(t(i), y(i));
        k2 = h * f(t(i) + h/2, y(i) + k1/2);
        k3 = h * f(t(i) + h/2, y(i) + k2/2);
        k4 = h * f(t(i) + h, y(i) + k3);

        y(i+1) = y(i) + (k1 + 2*k2 + 2*k3 + k4)/6;
    end
end

function y_prime = f(t, y)
    y_prime = t * y + t^3;
end

function y_exact = analyticalSolution(t)
    y_exact = exp(t^2 / 2 - 1);
end

rungeKuttaExample([0 1], 1, [5 10 20 40 80 160 320 640]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
