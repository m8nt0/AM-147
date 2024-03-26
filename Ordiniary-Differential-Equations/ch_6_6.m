function [t, y] = backwardEulerStiffExample()

    % Define the differential equation and initial condition
    f = @(t, y) y + 8*y^2 - 9*y^3;
    y0 = 1/2;

    % Define the time span
    tspan = [0, 3];

    % Set the step size
    h = 0.15;

    % Initialize arrays to store results
    t = tspan(1):h:tspan(2);
    y = zeros(size(t));

    % Set the initial condition
    y(1) = y0;

    % Solve using Backward Euler
    for i = 1:length(t)-1
        % Use Newton's method to solve for the next value
        y(i+1) = backwardEulerStep(t(i), y(i), h, f);
    end

    % Plot the results
    figure;
    plot(t, y, 'o-', 'DisplayName', 'Backward Euler');
    hold on;

    % Plot the true solution
    trueSolution = @(t) 1 ./ (1 + exp(-t));
    trueValues = trueSolution(t);
    plot(t, trueValues, '--', 'DisplayName', 'True Solution');

    xlabel('Time');
    ylabel('y(t)');
    title('Numerical Solution using Backward Euler');
    legend();

end

function y_new = backwardEulerStep(t, y, h, f)
    % Use Newton's method to solve for the next value in Backward Euler
    tol = 1e-8;
    maxIter = 100;
    z = y; % Initial guess

    for iter = 1:maxIter
        g = z - y - h * f(t + h, z);
        g_prime = 1 - h * (1 + 16*z - 27*z^2);
        z_new = z - g / g_prime;

        % Check for convergence
        if abs(z_new - z) < tol
            y_new = z_new;
            return;
        end

        z = z_new;
    end

    error('Newton''s method did not converge.');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Example usage
f = @(t, y) y + 8 * y^2 - 9 * y^3; % Your ODE function here
tspan = [0 3];
y0 = 1/2;
h = 0.15;
tol = 1e-8;

% [t, y] = backward_euler_solver(f, tspan, y0, h, tol);
% function [t, y] = backward_euler_solver(f, tspan, y0, h, tol)
t = tspan(1):h:tspan(2);
y = zeros(length(t), length(y0));
y(1, :) = y0;

for i = 1:length(t) - 1
    % Newton's Method to solve the implicit equation
    z = y(i, :); % Initial guess
    f_value = @(z) z - y(i, :) - h * (z + 8 * z^2 - 9 * z^3);
    df_value = @(z) 1 - h * (1 + 16 * z - 27 * z^2) - h * (16 * z + 8);
    znew = z - f_value(z) / df_value(z);
    
    % Continue Newton's Method until convergence
    while norm(znew - z) > tol
        z = znew;
        znew = z - f_value(z) / df_value(z);
    end
    
    % Update the solution
    y(i + 1, :) = znew;
end
% end


% Plot the solution
plot(t, y, '-o');
xlabel('t');
ylabel('y');
title('Backward Euler Method with Newton''s Method');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%