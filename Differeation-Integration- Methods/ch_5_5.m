% Additional Example 1: Gaussian Quadrature for ln(x)
% Approximate the integral of ln(x) from 1 to 2 using Gaussian Quadrature.

% Function to integrate
f = @(x) log(x);

% Transformation for [1,2] to [-1,1]
a = 1;
b = 2;
g = @(t) f((b - a) * t / 2 + (a + b) / 2);

% Degree of Gaussian Quadrature
n = 4;

% Nodes and weights for Gaussian Quadrature on [-1,1]
roots = [-0.86113631159405, -0.33998104358486, 0.33998104358486, 0.86113631159405];
weights = [0.34785484513745, 0.65214515486255, 0.65214515486255, 0.34785484513745];

% Gaussian Quadrature approximation
result = (b - a) / 2 * sum(weights .* g(roots));

% Display the result
disp(result);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Additional Example 2: Gaussian Quadrature for π
% Approximate π using the integral of 1/(1 + x^2) from 0 to 1 with Gaussian Quadrature.

% Function to integrate
f = @(x) 1 / (1 + x^2);

% Transformation for [0,1] to [-1,1]
a = 0;
b = 1;
g = @(t) f((b - a) * t / 2 + (a + b) / 2);

% Degree of Gaussian Quadrature
n = 4;

% Nodes and weights for Gaussian Quadrature on [-1,1]
roots = [-0.86113631159405, -0.33998104358486, 0, 0.33998104358486, 0.86113631159405];
weights = [0.34785484513745, 0.65214515486255, 0.88888888888888, 0.65214515486255, 0.34785484513745];

% Gaussian Quadrature approximation
result = (b - a) / 2 * sum(weights .* g(roots));

% Display the result
disp(result);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Function to represent the parametric path P = {x(t), y(t)}
x = @(t) 0.5 + 0.3*t + 3.9*t.^2 - 4.7*t.^3;
y = @(t) 1.5 + 0.3*t + 0.9*t.^2 - 2.7*t.^3;

% Define the derivative functions dx/dt and dy/dt
dx_dt = @(t) gradient(x(t), t);
dy_dt = @(t) gradient(y(t), t);

% Function for the arc length integrand
arc_length_integrand = @(t) sqrt(dx_dt(t).^2 + dy_dt(t).^2);

% Adaptive Quadrature function to compute arc length between t1 and t2
adap_quad_arc_length = @(t1, t2, tol) adapquad(arc_length_integrand, t1, t2, tol);

% Number of equal lengths (subpaths)
n = 10;

% Initialize variables to store results
arc_lengths = zeros(1, n);
t_values = zeros(1, n + 1);

% Calculate arc lengths and t values for equal intervals
for i = 1:n
    t1 = (i - 1) / n;
    t2 = i / n;
    
    % Apply Adaptive Quadrature
    arc_lengths(i) = adap_quad_arc_length(t1, t2, 1e-6);
    
    % Store t values
    t_values(i) = t1;
end
t_values(n + 1) = 1;  % Include the endpoint

% Display the results
disp('Equal Length Subpaths:');
disp('t-values:');
disp(t_values);
disp('Arc Lengths:');
disp(arc_lengths);
