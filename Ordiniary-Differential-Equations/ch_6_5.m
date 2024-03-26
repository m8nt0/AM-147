% Embedded Runge–Kutta Pairs
function [tout, yout, te, ye, ie] = embeddedRKPair(odefun, tspan, y0, options)
    % ode45 with embedded Runge–Kutta pairs
    
    % Set default options if not provided
    if nargin < 4
        options = odeset();
    end
    
    % Initialize variables
    tout = [];
    yout = [];
    te = [];
    ye = [];
    ie = [];
    
    % Event handling
    options.Events = @(t, y) events(t, y, odefun);
    
    % Initial step size
    h = (tspan(2) - tspan(1)) / 100;
    
    % Initialize
    t = tspan(1);
    y = y0(:);
    tout = t;
    yout = y.';
    
    while t < tspan(2)
        % Adjust the step size
        hmin = options.MinStep;
        hmax = options.MaxStep;
        h = min(hmax, max(hmin, h));
        
        % Take a step using an embedded Runge–Kutta pair
        [t, y, te_, ye_, ie_] = ode45(odefun, [t, t + h], y, options);
        
        % Append results
        tout = [tout; t(2:end)];
        yout = [yout; y(2:end, :)];
        
        % Store event information
        te = [te; te_];
        ye = [ye; ye_];
        ie = [ie; ie_];
        
        % Check for terminal events
        if ~isempty(ie_)
            break;
        end
        
        % Update the time and state for the next iteration
        t = t(end);
        y = y(end, :);
        
        % Adjust the step size based on error estimate
        h = adaptStepSize(options, h, y, y(end, :));
    end
end

% Event function
function [value, isterminal, direction] = events(t, y, odefun)
    value = odefun(t, y);
    isterminal = ones(size(value));
    direction = zeros(size(value));
end

% Adapt step size based on error estimate
function h = adaptStepSize(options, h, y, y_new)
    % Calculate error estimate
    err = max(abs(y_new - y) ./ max(abs(y), options.RelTol * options.AbsTol));
    
    % Adjust step size using a simple heuristic
    factor = 0.9 * (options.RelTol / err)^(1/5);
    
    % Limit the factor to avoid drastic changes
    factor = min(max(factor, 0.2), 5.0);
    
    % Update step size
    h = h * factor;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

odefun = @(t, y) -y;
tspan = [0, 5];
y0 = 1;
options = odeset('RelTol', 1e-5, 'MaxStep', 0.1);
[tout, yout] = embeddedRKPair(odefun, tspan, y0, options);

% Plot the results
figure;
plot(tout, yout);
xlabel('Time');
ylabel('y(t)');
title('ODE Solution using Embedded Runge–Kutta Pair');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [t, y] = bogacki_shampine_solver(f, tspan, y0, hmax, tol)
    t = tspan(1):hmax:tspan(2);
    y = zeros(length(t), length(y0));
    y(1, :) = y0;
    
    for i = 1:length(t) - 1
        h = min(hmax, tspan(2) - t(i));
        
        % Compute stages
        s1 = f(t(i), y(i, :));
        s2 = f(t(i) + 1/2 * h, y(i, :) + 1/2 * h * s1);
        s3 = f(t(i) + 3/4 * h, y(i, :) + 3/4 * h * s2);
        
        % Compute solutions
        zi1 = y(i, :) + h * (2/9 * s1 + 1/3 * s2 + 4/9 * s3);
        wi1 = y(i, :) + h * (7/24 * s1 + 1/4 * s2 + 1/3 * s3 + 1/8 * f(t(i) + h, zi1));
        
        % Compute error estimate
        ei1 = h * norm(-5/72 * s1 + 1/12 * s2 + 1/9 * s3 - 1/8 * f(t(i) + h, zi1));
        
        % Check if error is within tolerance
        if ei1 <= tol
            % Accept the step
            t(i + 1) = t(i) + h;
            y(i + 1, :) = wi1;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example usage
f = @(t, y) -y + t; % Your ODE function here
tspan = [0 2];
y0 = 1;
hmax = 0.1;
tol = 1e-5;

[t, y] = bogacki_shampine_solver(f, tspan, y0, hmax, tol);

% Plot the solution
plot(t, y, '-o');
xlabel('t');
ylabel('y');
title('Bogacki–Shampine ODE Solver');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [t, y] = runge_kutta_fehlberg_solver(f, tspan, y0, hmax, tol)
    t = tspan(1):hmax:tspan(2);
    y = zeros(length(t), length(y0));
    y(1, :) = y0;
    
    for i = 1:length(t) - 1
        h = min(hmax, tspan(2) - t(i));
        
        % Compute stages
        s1 = f(t(i), y(i, :));
        s2 = f(t(i) + 1/4 * h, y(i, :) + 1/4 * h * s1);
        s3 = f(t(i) + 3/10 * h, y(i, :) + 3/32 * h * s1 + 9/32 * h * s2);
        s4 = f(t(i) + 3/5 * h, y(i, :) + 1932/2197 * h * s1 - 7200/2197 * h * s2 + 7296/2197 * h * s3);
        s5 = f(t(i) + h, y(i, :) + 439/216 * h * s1 - 8 * h * s2 + 3680/513 * h * s3 - 845/4104 * h * s4);
        s6 = f(t(i) + 1/2 * h, y(i, :) - 8/27 * h * s1 + 2 * h * s2 - 3544/2565 * h * s3 + 1859/4104 * h * s4 - 11/40 * h * s5);
        
        % Compute solutions
        zi1 = y(i, :) + h * (35/384 * s1 + 500/1113 * s3 + 125/192 * s4 - 2187/6784 * s5 + 11/84 * s6);
        wi1 = y(i, :) + h * (5179/57600 * s1 + 7571/16695 * s3 + 393/640 * s4 - 92097/339200 * s5 + 187/2100 * s6 + 1/40 * f(t(i) + h, zi1));
        
        % Compute error estimate
        ei1 = h * norm(71/57600 * s1 - 71/16695 * s3 + 71/1920 * s4 - 17253/339200 * s5 + 22/525 * s6);
        
        % Check if error is within tolerance
        if ei1 <= tol
            % Accept the step
            t(i + 1) = t(i) + h;
            y(i + 1, :) = zi1;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example usage
f = @(t, y) 10 * (1 - y); % Your ODE function here
tspan = [0 100];
y0 = 0.5;
hmax = 1;
tol = 1e-4;

[t, y] = runge_kutta_fehlberg_solver(f, tspan, y0, hmax, tol);

% Plot the solution
plot(t, y, '-o');
xlabel('t');
ylabel('y');
title('Runge–Kutta–Fehlberg ODE Solver');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%