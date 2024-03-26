% Define the ODE function
ode_fun = @(t, y) -y + t^2 + 1;

% Time span
tspan = [0, 2];

% Initial condition
y0 = 0;

% Analyze Local and Global Truncation Error
h_values = [0.1, 0.05, 0.025];
for h = h_values
    % Explicit Euler Method (First-Order)
    [t_euler, y_euler] = euler_method(ode_fun, tspan, y0, h);

    % Global Truncation Error Analysis
    global_error = abs(y_euler(end) - exact_solution(t_euler(end)));
    fprintf('Global Truncation Error (h = %.3f): %.6f\n', h, global_error);
    
    % Plot results
    figure;
    plot(t_euler, y_euler, '-o', 'DisplayName', 'Euler Method');
    hold on;
    t_exact = linspace(tspan(1), tspan(2), 100);
    plot(t_exact, exact_solution(t_exact), 'k-', 'DisplayName', 'Exact Solution');
    xlabel('Time (t)');
    ylabel('Solution (y)');
    title(['Explicit Euler Method with h = ' num2str(h)]);
    legend('Location', 'NorthWest');
    grid on;
    hold off;
end

% Explicit Trapezoid Method
[t_trapezoid, y_trapezoid] = explicit_trapezoid_method(ode_fun, tspan, y0);

% Plot results
figure;
plot(t_trapezoid, y_trapezoid, '-o', 'DisplayName', 'Explicit Trapezoid Method');
hold on;
plot(t_exact, exact_solution(t_exact), 'k-', 'DisplayName', 'Exact Solution');
xlabel('Time (t)');
ylabel('Solution (y)');
title('Explicit Trapezoid Method');
legend('Location', 'NorthWest');
grid on;
hold off;

% Taylor Methods
order = 2; % Specify the order of Taylor method
[t_taylor, y_taylor] = taylor_method(ode_fun, tspan, y0, order);

% Plot results
figure;
plot(t_taylor, y_taylor, '-o', 'DisplayName', ['Taylor Method (Order ' num2str(order) ')']);
hold on;
plot(t_exact, exact_solution(t_exact), 'k-', 'DisplayName', 'Exact Solution');
xlabel('Time (t)');
ylabel('Solution (y)');
title(['Taylor Method (Order ' num2str(order) ')']);
legend('Location', 'NorthWest');
grid on;
hold off;

% Function to calculate the exact solution for comparison
function y = exact_solution(t)
    y = (t + 1).^2 - exp(-t);
end

% Explicit Euler Method (First-Order)
function [t, y] = euler_method(ode_fun, tspan, y0, h)
    t = tspan(1):h:tspan(2);
    y = zeros(size(t));
    y(1) = y0;
    for i = 1:length(t)-1
        y(i+1) = y(i) + h * ode_fun(t(i), y(i));
    end
end

% Explicit Trapezoid Method
function [t, y] = explicit_trapezoid_method(ode_fun, tspan, y0)
    t = tspan(1):0.1:tspan(2); % Use smaller step for plotting
    h = t(2) - t(1);
    y = zeros(size(t));
    y(1) = y0;
    for i = 1:length(t)-1
        f1 = ode_fun(t(i), y(i));
        f2 = ode_fun(t(i+1), y(i) + h * f1);
        y(i+1) = y(i) + h * (f1 + f2) / 2;
    end
end

% Taylor Method
function [t, y] = taylor_method(ode_fun, tspan, y0, order)
    t = tspan(1):0.1:tspan(2); % Use smaller step for plotting
    h = t(2) - t(1);
    y = zeros(size(t));
    y(1) = y0;
    for i = 1:length(t)-1
        f = zeros(1, order+1);
        for j = 0:order
            f(j+1) = ode_fun(t(i), y(i));
            y_deriv = zeros(1, j+1);
            y_deriv(end) = 1;
            for k = 1:j
                y_deriv(k) = (h^k) / factorial(k) * ode_fun(t(i), y(i));
            end
            f(j+1) = f(j+1) / factorial(j) * dot(y_deriv, f);
        end
        y(i+1) = y(i) + sum(f);
    end
end

%%%%%%%%%%%%%%  % %% % %%%%%%%

% Program 6.5: Euler’s Method and Global Truncation Error Calculation
% Input: interval inter, initial value y0, step size h
% Output: time steps t, Euler's Method approximate solution w, and global truncation error g
% Example usage: euler_global_error([-10 0], 1/10001, 1e-5);

function [t, w, g] = euler_global_error(inter, y0, h)
    % Define the function f(t, y) from the differential equation
    f = @(t, y) -4 * t^3 * y^2;

    % Initialize variables
    t = inter(1):h:inter(2);
    w = zeros(size(t));
    g = zeros(size(t));
    w(1) = y0;
    g(1) = 0; % Initial global truncation error is zero

    % Perform Euler's Method and calculate global truncation error
    for i = 1:length(t)-1
        w(i+1) = w(i) + h * f(t(i), w(i));
        z = 1 / (t(i+1)^4 + 1); % Exact solution at the next time step
        g(i+1) = abs(w(i+1) - z) + g(i); % Update global truncation error
    end

    % Plot the solution and global truncation error
    figure;
    subplot(2, 1, 1);
    plot(t, w);
    xlabel('t');
    ylabel('Approximate Solution w');
    title('Euler''s Method Approximation for Example 6.9');

    subplot(2, 1, 2);
    semilogy(t, g);
    xlabel('t');
    ylabel('Global Truncation Error g');
    title('Global Truncation Error for Example 6.9');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Program 6.6: Explicit Trapezoid Method
% Input: interval inter, initial value y0, step size h
% Output: time steps t and Trapezoid Method approximate solution w
% Example usage: trapezoid_method([-10 0], 1/10001, 1e-5);

function [t, w] = trapezoid_method(inter, y0, h)
    % Define the function f(t, y) from the differential equation
    f = @(t, y) -4 * t^3 * y^2;

    % Initialize variables
    t = inter(1):h:inter(2);
    w = zeros(size(t));
    w(1) = y0;

    % Perform Trapezoid Method
    for i = 1:length(t)-1
        % Calculate slopes
        sL = f(t(i), w(i));
        sR = f(t(i) + h, w(i) + h * sL);

        % Update solution using Trapezoid Method formula
        w(i+1) = w(i) + h/2 * (sL + sR);
    end

    % Plot the Trapezoid Method approximation
    figure;
    plot(t, w);
    xlabel('t');
    ylabel('Approximate Solution w');
    title('Explicit Trapezoid Method Approximation for Example 6.11');
end

trapezoid_method([-10 0], 1/10001, 1e-5);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Program for Additional Example 1: Trapezoid Method
% Input: interval inter, initial value y0, step size h
% Output: time steps t and Trapezoid Method approximate solution w

function [t, w] = trapezoid_method_additional_example1(inter, y0, h)
    % Define the function f(t, y) from the differential equation
    f = @(t, y) t * y^2;

    % Define the exact solution for error calculation
    exact_solution = @(t) -2 / (t^2 + 2);

    % Initialize variables
    t = inter(1):h:inter(2);
    w = zeros(size(t));
    w(1) = y0;

    % Perform Trapezoid Method
    for i = 1:length(t)-1
        % Calculate slopes
        sL = f(t(i), w(i));
        sR = f(t(i) + h, w(i) + h * sL);

        % Update solution using Trapezoid Method formula
        w(i+1) = w(i) + h/2 * (sL + sR);
    end

    % Plot the Trapezoid Method approximation
    figure;
    plot(t, w, '-o', 'DisplayName', 'Trapezoid Method');
    hold on;
    
    % Plot the exact solution
    exact_values = exact_solution(t);
    plot(t, exact_values, '-s', 'DisplayName', 'Exact Solution');
    
    xlabel('t');
    ylabel('Approximate Solution w / Exact Solution');
    title('Trapezoid Method Approximation for Additional Example 1');
    legend('show');
end

trapezoid_method_additional_example1([0 1], -1, 1/4);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

