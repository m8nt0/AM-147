% Chapter 5: Numerical Differentiation

% 5.1.1 Finite difference formulas

function derivative = numericalDifferentiation(func, x, method, h)
    % func: Function to be differentiated
    % x: Point at which differentiation is performed
    % method: 'forward', 'backward', or 'central'
    % h: Step size
    
    if nargin < 4
        h = 1e-5; % Default step size
    end
    
    switch method
        case 'forward'
            derivative = (func(x + h) - func(x)) / h;
        case 'backward'
            derivative = (func(x) - func(x - h)) / h;
        case 'central'
            derivative = (func(x + h) - func(x - h)) / (2 * h);
        otherwise
            error('Invalid differentiation method. Use ''forward'', ''backward'', or ''central''.');
    end
end

% Example of usage:
% Define your function
f = @(x) x^2;

% Choose a point
x0 = 2;

% Choose a differentiation method and compute the derivative
derivative_forward = numericalDifferentiation(f, x0, 'forward');
derivative_backward = numericalDifferentiation(f, x0, 'backward');
derivative_central = numericalDifferentiation(f, x0, 'central');

% Display results
disp(['Forward Difference: ' num2str(derivative_forward)]);
disp(['Backward Difference: ' num2str(derivative_backward)]);
disp(['Central Difference: ' num2str(derivative_central)]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Chapter 5: Numerical Differentiation

% 5.1.2 Rounding error

% Example to demonstrate rounding error
format long; % Display more decimal places

x = 0.1;
derivative_exact = -sin(x); % Exact derivative of sin(x) at x=0.1

h_values = 10.^(-1:-1:-15); % Various step sizes

for i = 1:length(h_values)
    derivative_approx = (sin(x + h_values(i)) - sin(x)) / h_values(i);
    rounding_error = abs(derivative_exact - derivative_approx);
    disp(['Step Size: ' num2str(h_values(i)) ', Rounding Error: ' num2str(rounding_error)]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Chapter 5: Numerical Differentiation

% 5.1.3 Extrapolation

% Example to demonstrate extrapolation
format long; % Display more decimal places

f = @(x) exp(x); % Example function
x0 = 0.1; % Point of interest
h = 0.1; % Step size

derivative_approx_h = (f(x0 + h) - f(x0)) / h;
derivative_approx_h_over_2 = (f(x0 + h/2) - f(x0 - h/2)) / (h/2);

extrapolated_derivative = (4 * derivative_approx_h_over_2 - derivative_approx_h) / 3;

disp(['Exact Derivative: ' num2str(f(x0))]);
disp(['Extrapolated Derivative: ' num2str(extrapolated_derivative)]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Chapter 5: Numerical Differentiation

% 5.1.4 Symbolic differentiation and integration

syms x;

% Symbolic differentiation
f = x^2;
derivative_symbolic = diff(f, x);
disp(['Symbolic Derivative: ' char(derivative_symbolic)]);

% Symbolic integration
integral_symbolic = int(f, x);
disp(['Symbolic Integral: ' char(integral_symbolic)]);
