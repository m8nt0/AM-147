function [t, y] = euler_method(f, tspan, y0, h)
    % Euler's method for solving an ODE: dy/dt = f(t, y)
    % Input:
    %   - f: function handle for the ODE, f(t, y)
    %   - tspan: time span [t_start, t_end]
    %   - y0: initial condition at t_start
    %   - h: step size
    % Output:
    %   - t: time values
    %   - y: corresponding solution values

    % Initialize arrays
    t = tspan(1):h:tspan(2);
    y = zeros(size(t));

    % Set initial condition
    y(1) = y0;

    % Euler's method iteration
    for i = 1:length(t)-1
        y(i+1) = y(i) + h * f(t(i), y(i));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define coefficient functions a(t) and b(t)
a = @(t) -0.5; % You can change this function as needed
b = @(t) 2 * sin(t); % You can change this function as needed

% Time span
tspan = [0, 5];

% Initial condition
y0 = 1;

% Define the ODE function
ode_fun = @(t, y) a(t) * y + b(t);

% Solve the ODE using ode45 (MATLAB's built-in solver)
[t, y] = ode45(ode_fun, tspan, y0);

% Plot the solution
figure;
plot(t, y, 'b-', 'LineWidth', 2);
xlabel('Time (t)');
ylabel('Solution (y)');
title('Solution to \frac{dy}{dt} = a(t)y + b(t)');
grid on;

% Display a message about existence and uniqueness
disp('Existence and Uniqueness:');
disp('For this example, the ODE has a continuous and bounded solution for any t.');
disp('The coefficients a(t) and b(t) are assumed to satisfy the conditions for existence and uniqueness.');

% Check continuity of the solution
cont_check = iscontinuous(t, y);
disp(['Continuity Check: The solution is continuous on the given time span? ' num2str(cont_check)]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [t, y] = solve_linear_ode(a, b, tspan, y0)
    % Solve a first-order linear ODE: dy/dt = a(t)y + b(t)
    % Input:
    %   - a: function handle for coefficient a(t)
    %   - b: function handle for coefficient b(t)
    %   - tspan: time span [t_start, t_end]
    %   - y0: initial condition at t_start
    % Output:
    %   - t: time values
    %   - y: corresponding solution values

    % Define the ODE function
    f = @(t, y) a(t) * y + b(t);

    % Solve using ode45 (built-in MATLAB solver)
    [t, y] = ode45(f, tspan, y0);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Program 6.1: Euler’s Method for Solving Initial Value Problems
% Use with ydot.m to evaluate rhs of differential equation
% Input: interval inter, initial value y0, number of steps n
% Output: time steps t, solution y
% Example usage: euler1([0 1], 1, 10);

function [t, y] = euler1(inter, y0, n)
    t(1) = inter(1);
    y(1) = y0;
    h = (inter(2) - inter(1)) / n;

    for i = 1:n
        t(i+1) = t(i) + h;
        y(i+1) = eulerstep(t(i), y(i), h);
    end

    % Plotting the solution
    plot(t, y)
    xlabel('t')
    ylabel('y')
    title('Euler''s Method for y'' = ty + t^3')
end

function y = eulerstep(t, y, h)
    % One step of Euler’s Method
    % Input: current time t, current value y, stepsize h
    % Output: approximate solution value at time t+h
    y = y + h * ydot(t, y);
end

function z = ydot(t, y)
    % Right-hand side of differential equation
    z = t * y + t^3;
end

euler1([0 1], 1, 10);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Program 6.2: Euler’s Method for Initial Value Problems with Error Estimation
% Input: interval inter, initial value y0, number of steps n
% Output: time steps t, solution y, estimated error est
% Example usage: euler2([0 1], 1, 10);

function [t, y, est] = euler2(inter, y0, n)
    t(1) = inter(1);
    y(1) = y0;
    h = (inter(2) - inter(1)) / n;

    for i = 1:n
        t(i+1) = t(i) + h;
        [y(i+1), est(i)] = eulerstep(t(i), y(i), h);
    end

    % Plotting the solution and estimated error
    figure;
    subplot(2,1,1);
    plot(t, y);
    xlabel('t');
    ylabel('y');
    title('Euler''s Method for y'' = ty + t^3');

    subplot(2,1,2);
    plot(t, est);
    xlabel('t');
    ylabel('Estimated Error');
    title('Estimated Error in Euler''s Method');
end

function [y, est] = eulerstep(t, y, h)
    % One step of Euler’s Method with error estimation
    % Input: current time t, current value y, stepsize h
    % Output: approximate solution value at time t+h, estimated error
    y_new = y + h * ydot(t, y);
    est = h^2 * ydot2(t, y_new) / 2;
    y = y_new;
end

function z = ydot(t, y)
    % Right-hand side of differential equation
    z = t * y + t^3;
end

function z = ydot2(t, y)
    % Second derivative of y with respect to t
    z = t + 3 * t^2;
end

euler2([0 1], 1, 10);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Program 6.3: Solving First-Order Linear Differential Equation
% Input: interval inter, initial value y0
% Output: time steps t, solution y
% Example usage: solve_linear_eqn([0 1], 1);

function [t, y] = solve_linear_eqn(inter, y0)
    % Define the function g(t) and h(t) from the differential equation
    g = @(t) t;
    h = @(t) t^3;

    % Calculate the integrating factor
    int_factor = exp(-integral(g, inter(1), inter(2)));

    % Define the differential equation
    dydt = @(t, y) g(t) * y + h(t);

    % Use the integrating factor to find the solution
    y_sol = @(t) int_factor * integral(@(u) exp(integral(g, inter(1), u)) * h(u), inter(1), t);

    % Generate time steps and corresponding solution values
    t = linspace(inter(1), inter(2), 100);
    y = arrayfun(y_sol, t);

    % Plot the solution
    figure;
    plot(t, y);
    xlabel('t');
    ylabel('y');
    title('Solution of First-Order Linear Differential Equation');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Program 6.4: Euler’s Method for First-Order Linear Differential Equation
% Input: interval inter, initial value y0, step size h
% Output: time steps t, Euler's Method approximate solution y
% Example usage: euler_linear_eqn([0 1], 1, 0.1);

function [t, y] = euler_linear_eqn(inter, y0, h)
    t(1) = inter(1);
    y(1) = y0;

    % Define the function g(t) and h(t) from the differential equation
    g = @(t) t;
    h_t = @(t) t^3;

    % Perform Euler's Method
    for i = 1:((inter(2) - inter(1)) / h)
        t(i+1) = t(i) + h;
        y(i+1) = y(i) + h * (g(t(i)) * y(i) + h_t(t(i)));
    end

    % Plot the Euler's Method approximate solution
    figure;
    plot(t, y);
    xlabel('t');
    ylabel('Euler''s Method Approximation');
    title('Euler''s Method Approximation for First-Order Linear Differential Equation');
end


solve_linear_eqn([0 1], 1);

euler_linear_eqn([0 1], 1, 0.1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

