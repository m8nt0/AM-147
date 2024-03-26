% Runge-Kutta Methods

% Define the ODE system
ode_fun = @(t, Y) -2 * Y;

% Time span
tspan = [0, 5];

% Initial condition
initial_condition = 1;

% Solve using various Runge-Kutta methods
[t_rk1, Y_rk1] = ode45(ode_fun, tspan, initial_condition); % RK45
[t_rk2, Y_rk2] = ode23(ode_fun, tspan, initial_condition); % RK23
[t_rk4, Y_rk4] = ode113(ode_fun, tspan, initial_condition); % RK113

% Plot the results
figure;
plot(t_rk1, Y_rk1, '-o', 'DisplayName', 'RK45');
hold on;
plot(t_rk2, Y_rk2, '-o', 'DisplayName', 'RK23');
plot(t_rk4, Y_rk4, '-o', 'DisplayName', 'RK113');
xlabel('Time (t)');
ylabel('Solution');
title('Runge-Kutta Methods Comparison');
legend;
grid on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Hodgkin-Huxley Neuron Simulation

% Define the ODE system for Hodgkin-Huxley model
hh_ode = @(t, Y) hh_neuron_ode(t, Y);

% Time span
tspan_hh = [0, 20];

% Initial conditions for the Hodgkin-Huxley model
initial_conditions_hh = hh_neuron_initial_conditions();

% Solve the system
[t_hh, Y_hh] = ode45(hh_ode, tspan_hh, initial_conditions_hh);

% Plot the results
figure;
plot(t_hh, Y_hh(:, 1), '-o', 'DisplayName', 'Voltage (V)');
xlabel('Time (ms)');
ylabel('Voltage (mV)');
title('Hodgkin-Huxley Neuron Simulation');
legend;
grid on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Lorenz Equations Simulation

% Define the ODE system for Lorenz equations
lorenz_ode = @(t, Y) lorenz_equations(t, Y);

% Time span
tspan_lorenz = [0, 25];

% Initial conditions for the Lorenz equations
initial_conditions_lorenz = [0, 1, 0];

% Solve the system
[t_lorenz, Y_lorenz] = ode45(lorenz_ode, tspan_lorenz, initial_conditions_lorenz);

% Plot the results
figure;
plot3(Y_lorenz(:, 1), Y_lorenz(:, 2), Y_lorenz(:, 3), 'LineWidth', 1.5);
xlabel('X-axis');
ylabel('Y-axis');
zlabel('Z-axis');
title('Lorenz Equations Simulation');
grid on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Program 6.6 Animation program for bridge using IVP solver
% Inputs: inter = time interval inter,
% ic = [y(1,1) y(1,2) y(1,3) y(1,4)],
% number of steps n, steps per point plotted p
% Calls a one-step method such as trapstep.m
% Example usage: tacoma([0 1000],[1 0 0.001 0],25000,5);

function tacoma(inter, ic, n, p)
    clf % clear figure window
    h = (inter(2) - inter(1)) / n;
    y(1, :) = ic; % enter initial conds in y
    t(1) = inter(1); len = 6;
    set(gca, 'XLim', [-8 8], 'YLim', [-8 8], ...
        'XTick', [-8 0 8], 'YTick', [-8 0 8]);
    cla; % clear screen
    axis square % make aspect ratio 1 - 1
    road = animatedline('color', 'b', 'LineStyle', '-', 'LineWidth', 1);
    lcable = animatedline('color', 'r', 'LineStyle', '-', 'LineWidth', 1);
    rcable = animatedline('color', 'r', 'LineStyle', '-', 'LineWidth', 1);
    
    for k = 1:n
        for i = 1:p
            t(i+1) = t(i) + h;
            y(i+1, :) = trapstep(t(i), y(i, :), h);
        end
        y(1, :) = y(p+1, :);
        t(1) = t(p+1);
        z1(k) = y(1, 1);
        z3(k) = y(1, 3);
        c = len * cos(y(1, 3));
        s = len * sin(y(1, 3));
        clearpoints(road);
        addpoints(road, [-c c], [-s-y(1, 1) s-y(1, 1)])
        clearpoints(lcable);
        addpoints(lcable, [-c -c], [-s-y(1, 1) 8])
        clearpoints(rcable);
        addpoints(rcable, [c c], [s-y(1, 1) 8])
        drawnow;
        pause(h)
    end
end

function y = trapstep(t, x, h)
    % one step of the Trapezoid Method
    z1 = ydot(t, x);
    g = x + h * z1;
    z2 = ydot(t + h, g);
    y = x + h * (z1 + z2) / 2;
end

function ydot = ydot(t, y)
    len = 6; a = 0.2; W = 80; omega = 2 * pi * 38 / 60;
    a1 = exp(a * (y(1) - len * sin(y(3))));
    a2 = exp(a * (y(1) + len * sin(y(3))));
    ydot(1) = y(2);
    ydot(2) = -0.01 * y(2) - 0.4 * (a1 + a2 - 2) / a + 0.2 * W * sin(omega * t);
    ydot(3) = y(4);
    ydot(4) = -0.01 * y(4) + 1.2 * cos(y(3)) * (a1 - a2) / (len * a);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Program for Example 6.5: Hodgkin-Huxley neuron simulation
% Example usage: hh_simulation([0,100],[-65,0,0.3,0.6],2000);

function hh_simulation(inter, ic, n)
    global pa pb pulse
    fprintf('Enter pulse start, end, muamps in [ ], e.g. [50 51 7]: ');
    inp = input('');
    pa = inp(1);
    pb = inp(2);
    pulse = inp(3);
    
    a = inter(1);
    b = inter(2);
    h = (b - a) / n; % plot n points in total
    y(1, :) = ic; % enter initial conds in y
    t(1) = a;
    
    for i = 1:n
        t(i+1) = t(i) + h;
        y(i+1, :) = rk4step(t(i), y(i, :), h);
    end
    
    subplot(3,1,1);
    plot([a pa pa pb pb b], [0 0 pulse pulse 0 0]);
    grid; axis([0 100 0 2*pulse]);
    ylabel('input pulse');
    
    subplot(3,1,2);
    plot(t, y(:, 1));
    grid; axis([0 100 -100 100]);
    ylabel('voltage (mV)');
    
    subplot(3,1,3);
    plot(t, y(:, 2), t, y(:, 3), t, y(:, 4));
    grid; axis([0 100 0 1]);
    ylabel('gating variables');
    legend('m', 'n', 'h');
    
    xlabel('time (msec)');
end

function y = rk4step(t, w, h)
    % one step of the Runge-Kutta order 4 method
    s1 = ydot(t, w);
    s2 = ydot(t + h/2, w + h*s1/2);
    s3 = ydot(t + h/2, w + h*s2/2);
    s4 = ydot(t + h, w + h*s3);
    y = w + h * (s1 + 2*s2 + 2*s3 + s4) / 6;
end

function z = ydot(t, w)
    global pa pb pulse
    c = 1; g1 = 120; g2 = 36; g3 = 0.3; T = (pa + pb) / 2; len = pb - pa;
    e0 = -65; e1 = 50; e2 = -77; e3 = -54.4;
    in = pulse * (1 - sign(abs(t - T) - len/2)) / 2;
    
    v = w(1); m = w(2); n = w(3); h = w(4);
    z = zeros(1, 4);
    z(1) = (in - g1*m^3*h*(v - e1) - g2*n^4*(v - e2) - g3*(v - e3)) / c;
    
    v = v - e0;
    z(2) = (1 - m) * (2.5 - 0.1*v) / (exp(2.5 - 0.1*v) - 1) - m * 4 * exp(-v/18);
    z(3) = (1 - n) * (0.1 - 0.01*v) / (exp(1 - 0.1*v) - 1) - n * 0.125 * exp(-v/80);
    z(4) = (1 - h) * 0.07 * exp(-v/20) - h / (exp(3 - 0.1*v) + 1);
end

hh_simulation([0 100], [-65 0 0.3 0.6], 2000);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Lorenz equations
function z = lorenz(t, y)
    s = 10; r = 28; b = 8/3;
    z(1) = -s * y(1) + s * y(2);
    z(2) = -y(1) * y(3) + r * y(1) - y(2);
    z(3) = y(1) * y(2) - b * y(3);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Example 1: y = t*y, y(0) = 1
h = 1/4;
tspan = [0, 1];
y0 = 1;
[t, w] = rungeKuttaOrder4(@example1ODE, tspan, y0, h);
globalTruncationError = abs(exp(1) - w(end));

% Display result
fprintf('Approximate solution at t=1: %f\n', w(end));
fprintf('Global Truncation Error at t=1: %f\n', globalTruncationError);

% Function for Example 1 ODE
function dydt = example1ODE(t, y)
    dydt = t * y;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Example 2: y = e^sin(t), y(0) = 1
h = 1/2;
tspan = [0, 4];
y0 = [exp(sin(0)); cos(0) * exp(sin(0))];
[t, w] = rungeKuttaOrder4(@example2ODE, tspan, y0, h);
globalTruncationError = abs(exp(sin(4)) - w(end, 1));

% Display result
figure;
plot(t, w(:, 1), '-o', 'DisplayName', 'Approximate Solution');
xlabel('t');
ylabel('y');
legend('Location', 'Best');
title('Runge-Kutta Order 4 Approximate Solution for Example 2');

% Function for Example 2 ODE
function dydt = example2ODE(t, y)
    dydt = [y(2); y(2) * sin(t) - y(1) * cos(t)];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Program 6.6 Animation program for bridge using IVP solver
% Inputs: inter = time interval inter,
% ic = [y(1,1) y(1,2) y(1,3) y(1,4)],
% number of steps n, steps per point plotted p
% Calls a one-step method such as trapstep.m
% Example usage: tacoma([0 1000],[1 0 0.001 0],25000,5);

function tacoma(inter, ic, n, p)
    clf % clear figure window
    h = (inter(2) - inter(1)) / n;
    y(1, :) = ic; % enter initial conds in y
    t(1) = inter(1);
    len = 6;
    set(gca, 'XLim', [-8 8], 'YLim', [-8 8], ...
        'XTick', [-8 0 8], 'YTick', [-8 0 8]);
    cla; % clear screen
    axis square % make aspect ratio 1 - 1
    road = animatedline('color', 'b', 'LineStyle', '-', 'LineWidth', 1);
    lcable = animatedline('color', 'r', 'LineStyle', '-', 'LineWidth', 1);
    rcable = animatedline('color', 'r', 'LineStyle', '-', 'LineWidth', 1);

    for k = 1:n
        for i = 1:p
            t(i + 1) = t(i) + h;
            y(i + 1, :) = trapstep(t(i), y(i, :), h);
        end
        y(1, :) = y(p + 1, :);
        t(1) = t(p + 1);
        z1(k) = y(1, 1);
        z3(k) = y(1, 3);
        c = len * cos(y(1, 3));
        s = len * sin(y(1, 3));
        clearpoints(road);
        addpoints(road, [-c c], [-s - y(1, 1) s - y(1, 1)]);
        clearpoints(lcable);
        addpoints(lcable, [-c -c], [-s - y(1, 1) 8]);
        clearpoints(rcable);
        addpoints(rcable, [c c], [s - y(1, 1) 8]);
        drawnow;
        pause(h)
    end

    function y = trapstep(t, x, h)
        %one step of the Trapezoid Method
        z1 = ydot(t, x);
        g = x + h * z1;
        z2 = ydot(t + h, g);
        y = x + h * (z1 + z2) / 2;
    end

    function ydot = ydot(t, y)
        len = 6;
        a = 0.2;
        W = 80;
        omega = 2 * pi * 38 / 60;
        a1 = exp(a * (y(1) - len * sin(y(3))));
        a2 = exp(a * (y(1) + len * sin(y(3))));
        ydot(1) = y(2);
        ydot(2) = -0.01 * y(2) - 0.4 * (a1 + a2 - 2) / a + 0.2 * W * sin(omega * t);
        ydot(3) = y(4);
        ydot(4) = -0.01 * y(4) + 1.2 * cos(y(3)) * (a1 - a2) / (len * a);
    end
end

tacoma([0 1000],[1 0 0.001 0],25000,5);
