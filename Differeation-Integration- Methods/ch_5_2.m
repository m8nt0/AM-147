% Chapter 5: Numerical Integration

% 5.2.1 Trapezoid Rule

% Example of Trapezoid Rule for numerical integration
f = @(x) x^2; % Example function to integrate
a = 0; % Lower limit of integration
b = 1; % Upper limit of integration
n = 4; % Number of subintervals

h = (b - a) / n; % Step size

x_values = a:h:b;
y_values = f(x_values);

integral_approx_trapezoid = h * (sum(y_values) - 0.5 * (y_values(1) + y_values(end)));

disp(['Approximate Integral using Trapezoid Rule: ' num2str(integral_approx_trapezoid)]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Chapter 5: Numerical Integration

% 5.2.2 Simpson’s Rule

% Example of Simpson's Rule for numerical integration
f = @(x) x^2; % Example function to integrate
a = 0; % Lower limit of integration
b = 1; % Upper limit of integration
n = 4; % Number of subintervals (must be even)

h = (b - a) / n; % Step size

x_values = a:h:b;
y_values = f(x_values);

integral_approx_simpson = h/3 * (y_values(1) + 4*sum(y_values(2:2:end-1)) + 2*sum(y_values(3:2:end-2)) + y_values(end));

disp(['Approximate Integral using Simpson’s Rule: ' num2str(integral_approx_simpson)]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Chapter 5: Numerical Integration

% 5.2.3 Composite Newton-Cotes formulas

% Example of Composite Trapezoid Rule for numerical integration
f = @(x) x^2; % Example function to integrate
a = 0; % Lower limit of integration
b = 1; % Upper limit of integration
n = 4; % Number of subintervals

h = (b - a) / n; % Step size

x_values = a:h:b;
y_values = f(x_values);

integral_approx_composite_trapezoid = h/2 * (y_values(1) + 2*sum(y_values(2:end-1)) + y_values(end));

disp(['Approximate Integral using Composite Trapezoid Rule: ' num2str(integral_approx_composite_trapezoid)]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Chapter 5: Numerical Integration

% 5.2.4 Open Newton–Cotes Methods

% Example of Open Trapezoid Rule for numerical integration
f = @(x) x^2; % Example function to integrate
a = 0; % Lower limit of integration
b = 1; % Upper limit of integration
n = 4; % Number of subintervals

h = (b - a) / (n + 2); % Step size

x_values = a + h:h:b - h;
y_values = f(x_values);

integral_approx_open_trapezoid = h * sum(y_values);

disp(['Approximate Integral using Open Trapezoid Rule: ' num2str(integral_approx_open_trapezoid)]);
