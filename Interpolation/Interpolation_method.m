function y = lagrange_interpolation(x_data, y_data, x)
    n = length(x_data);
    y = 0;

    for i = 1:n
        term = y_data(i);
        for j = 1:n
            if j ~= i
                term = term * (x - x_data(j)) / (x_data(i) - x_data(j));
            end
        end
        y = y + term;
    end
end

function y = newton_divided_differences(x_data, y_data, x)
    n = length(x_data);
    f = zeros(n, n);
    f(:, 1) = y_data';

    for j = 2:n
        for i = 1:n-j+1
            f(i, j) = (f(i+1, j-1) - f(i, j-1)) / (x_data(i+j-1) - x_data(i));
        end
    end

    y = 0;
    for j = 1:n
        term = f(1, j);
        for i = 1:j-1
            term = term * (x - x_data(i));
        end
        y = y + term;
    end
end

function result_lagrange_newton = interpolation_examples(x_data, y_data, x_interpolate)
    % Lagrange interpolation
    result_lagrange = lagrange_interpolation(x_data, y_data, x_interpolate);

    % Newton's divided differences
    result_newton = newton_divided_differences(x_data, y_data, x_interpolate);

    result_lagrange_newton = [result_lagrange, result_newton];
end

function P = lagrange_interpolation(x, y, xi)
    n = length(x);
    P = 0;

    for i = 1:n
        term = y(i);
        for j = 1:n
            if j ~= i
                term = term * (xi - x(j)) / (x(i) - x(j));
            end
        end
        P = P + term;
    end
end

function c = newton_divided_differences(x, y)
    n = length(x);
    v = zeros(n, n);

    for j = 1:n
        v(j, 1) = y(j);
    end

    for i = 2:n
        for j = 1:n + 1 - i
            v(j, i) = (v(j + 1, i - 1) - v(j, i - 1)) / (x(j + i - 1) - x(j));
        end
    end

    c = v(1, :);
end

function clickinterp
    xl = -3;
    xr = 3;
    yb = -3;
    yt = 3;
    plot([xl xr], [0 0], 'k', [0 0], [yb yt], 'k');
    grid on;
    
    xlist = [];
    ylist = [];
    k = 0;

    while (0 == 0)
        [xnew, ynew] = ginput(1);
        
        if length(xnew) < 1
            break;
        end
        
        k = k + 1;
        xlist(k) = xnew;
        ylist(k) = ynew;
        c = newtdd(xlist, ylist, k);
        x = xl:0.01:xr;
        y = nest(k - 1, c, x, xlist);
        plot(xlist, ylist, 'o', x, y, [xl xr], [0, 0], 'k', [0 0], [yb yt], 'k');
        axis([xl xr yb yt]);
        grid on;
    end
end

function y = sin_approximation(x)
    b = pi * (0:3) / 6;
    yb = sin(b);
    c = newtdd(b, yb, 4);

    s = 1;
    x1 = mod(x, 2 * pi);

    if x1 > pi
        x1 = 2 * pi - x1;
        s = -1;
    end

    if x1 > pi / 2
        x1 = pi - x1;
    end

    y = s * nest(3, c, x1, b);
end

function error = interpolation_error(x_data, y_data, x_interpolate, true_function)
    % Evaluate the true function at the interpolation point
    true_value = true_function(x_interpolate);
    
    % Perform interpolation
    interpolated_value = newton_divided_differences(x_data, y_data, x_interpolate);

    % Calculate interpolation error
    error = abs(interpolated_value - true_value);
end

syms x;
syms f(x);
n = 3; % Change n based on the degree of the polynomial

% Define a symbolic function f(x)
f(x) = x^2 + 2*x + 1;

% Generate symbolic variables for the interpolation points and coefficients
syms x_i;
syms a;

% Generate Newton's form
newton_form = a;
for i = 1:n
    newton_form = newton_form + a * prod(x_i - sym('x', [1:i-1])) * (x - x_i);
end

disp('Newton''s Form:');
disp(newton_form);

% Generate the error formula
error_formula = f(x) - newton_form;

disp('Error Formula:');
disp(error_formula);

% Simplify the expressions
newton_form_simplified = simplify(newton_form);
error_formula_simplified = simplify(error_formula);

disp('Newton''s Form (Simplified):');
disp(newton_form_simplified);

disp('Error Formula (Simplified):');
disp(error_formula_simplified);

function runge_phenomenon_example(n)
    % Function to demonstrate Runge phenomenon with n interpolation points
    x_data = linspace(-5, 5, n);
    y_data = 1./(1 + x_data.^2);

    % Interpolation points for plotting
    x_interpolate = linspace(-5, 5, 1000);

    % Interpolation using Newton's divided differences
    y_interpolate = zeros(size(x_interpolate));
    for i = 1:length(x_interpolate)
        y_interpolate(i) = newton_divided_differences(x_data, y_data, x_interpolate(i));
    end

    % Plotting
    figure;
    plot(x_interpolate, y_interpolate, 'r-', x_data, y_data, 'bo');
    title('Runge Phenomenon Example');
    xlabel('x');
    ylabel('y');
    legend('Interpolation', 'True Function');
end

function error_bound = calculate_interpolation_error_bound(x, xi)
    n = length(xi);
    product_term = prod(x - xi);
    f_nth_derivative = sin(x);  % Adjust this based on your actual function
    
    error_bound = abs(product_term / factorial(n) * f_nth_derivative);
end

% Runge Phenomenon Interpolation
function interpolated_values = interpolate_runge_function(xi_values)
    f = @(x) 1 ./ (1 + 12*x.^2);  % Runge function
    interpolated_values = polyval(polyfit(xi_values, f(xi_values), length(xi_values) - 1), xi_values);
end

% Chebyshev Interpolation
function y = chebyshev_interpolation(f, n, a, b, x)
    % Chebyshev Interpolation
    cheb_nodes = cos((2*(1:n)-1)*pi/(2*n));
    scaled_nodes = 0.5*(b-a)*cheb_nodes + 0.5*(b+a);
    f_values = f(scaled_nodes);
    
    % Coefficients using Chebyshev polynomials
    coeffs = chebyshev_coeffs(f_values);
    
    % Evaluate interpolation
    y = chebyshev_eval(coeffs, x, a, b);
end

function coeffs = chebyshev_coeffs(f_values)
    % Calculate Chebyshev coefficients
    n = length(f_values);
    coeffs = zeros(1, n);
    
    for k = 1:n
        coeffs(k) = (2/n) * sum(f_values .* cos((k-1)*acos(linspace(-1, 1, n))));
    end
end

function y = chebyshev_eval(coeffs, x, a, b)
    % Evaluate Chebyshev interpolation
    n = length(coeffs);
    y = zeros(size(x));
    
    for k = 1:n
        y = y + coeffs(k) * cos((k-1)*acos((2*(x-a)/(b-a))-1));
    end
end

function T = chebyshev_polynomials(n, x)
    % Generate Chebyshev polynomials up to degree n
    T = zeros(n+1, length(x));
    T(1, :) = ones(1, length(x));
    T(2, :) = x;
    
    for k = 3:n+1
        T(k, :) = 2 * x .* T(k-1, :) - T(k-2, :);
    end
end

function y_new = change_of_interval(y, a_old, b_old, a_new, b_new)
    % Change of interval for Chebyshev interpolation result
    y_new = ((b_new - a_new) / (b_old - a_old)) * (y - a_old) + a_new;
end

% Cubic Spline Interpolation
function y_interp = cubic_spline_interpolation(x_data, y_data, x_interp)
    n = length(x_data);
    h = diff(x_data);
    alpha = 3./h .* (diff(y_data)./h - diff(y_data([1 1:n-1]))./diff(x_data([1 1:n-1])));
    l = [1; 2*(h(1:n-2)+h(2:n-1)); 1];
    mu = [0; h(2:n-1); 0];
    z = tridiagonal_solver(l, diag([2; 2*(h(1:n-2)+h(2:n-1)); 2]), mu, alpha);

    a = y_data;
    b = (1./h).*(diff(y_data)./h - h.*(2*z(1:n-1)+z(2:n))/3);
    c = z(1:n-1);
    d = (z(2:n) - z(1:n-1))./(3*h);
    
    % Interpolation
    y_interp = zeros(size(x_interp));
    for i = 1:n-1
        indices = x_interp >= x_data(i) & x_interp <= x_data(i+1);
        x_interval = x_interp(indices) - x_data(i);
        y_interp(indices) = a(i) + b(i)*x_interval + c(i)*x_interval.^2 + d(i)*x_interval.^3;
    end
end

function y_interp = cubic_spline_endpoint_conditions(x_data, y_data, endpoint_condition)
    n = length(x_data);
    h = diff(x_data);
    alpha = 3./h .* (diff(y_data)./h - diff(y_data([1 1:n-1]))./diff(x_data([1 1:n-1])));
    
    if strcmp(endpoint_condition, 'natural')
        l = [1; 2*(h(1:n-2)+h(2:n-1)); 1];
        mu = [0; h(2:n-1); 0];
        z = tridiagonal_solver(l, diag([2; 2*(h(1:n-2)+h(2:n-1)); 2]), mu, alpha);
    elseif strcmp(endpoint_condition, 'clamped')
        h_0 = x_data(2) - x_data(1);
        h_n = x_data(n) - x_data(n-1);
        alpha_0 = (3/h_0)*(y_data(2) - y_data(1)) - (3/h)alpha(1);
        alpha_n = (3/h_n)(y_data(n) - y_data(n-1)) - (3/h)*alpha(n-1);
        l = [2*h_0 h(1:n-2)+h(2:n-1) 2*h_n];
        mu = [h_0 h(2:n-1) h_n];
        alpha = [alpha_0; alpha; alpha_n];
    
        z = tridiagonal_solver(l, diag([2*h_0 h(1:n-2)+2*(h(2:n-1)) 2*h_n]), mu, alpha);
    else
        error('Invalid endpoint condition. Choose ''natural'' or ''clamped''.');
    end
    
    a = y_data;
    b = (1./h).*(diff(y_data)./h - h.*(2*z(1:n-1)+z(2:n))/3);
    c = z(1:n-1);
    d = (z(2:n) - z(1:n-1))./(3*h);
    
    % Interpolation
    y_interp = zeros(size(x_data));
    for i = 1:n-1
        indices = x_interp >= x_data(i) & x_interp <= x_data(i+1);
        x_interval = x_interp(indices) - x_data(i);
        y_interp(indices) = a(i) + b(i)*x_interval + c(i)*x_interval.^2 + d(i)*x_interval.^3;
    end
end

function coeff = naturalCubicSplineCoefficients(x, y)
    n = length(x);
    A = zeros(n, n);
    r = zeros(n, 1);
    for i = 2:n-1
        dx = x(i) - x(i-1);
        dy = y(i) - y(i-1);
        A(i, i-1:i+1) = [dx, 2 * (dx + (x(i+1) - x(i))), x(i+1) - x(i)];
        r(i) = 3 * (dy / dx - (y(i) - y(i-1)) / (x(i) - x(i-1)));
    end
    
    % Natural spline conditions
    A(1, 1) = 2;
    A(n, n) = 2;
    r(1) = 0;
    r(n) = 0;
    
    % Solve the system of equations
    coeff = A \ r;
    
    % Adjust the number of coefficients
    coeff = [zeros(1, 3); coeff'; zeros(1, 3)];
end

x = [1, 2, 4, 5];
y = [2, 1, 4, 3];

coefficients = naturalCubicSplineCoefficients(x, y);
disp(coefficients);


